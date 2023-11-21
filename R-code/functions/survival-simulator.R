### survival-simulator.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Aug 19 2022 (13:51) 
## Version: 
## Last-Updated: Nov 12 2023 (14:59) 
##           By: Anders Munch
##     Update #: 92
#----------------------------------------------------------------------
## 
### Commentary: Flexible setup for simulating from general survival models
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(here)
library(data.table)
library(prodlim)
library(riskRegression)
library(ggplot2)
library(gridExtra)
## ----------------- Main function -----------------
## Main function to build survival simulator from smaller functions below
build_sim <- function(outcome, cens, covars = sim_covars_truncated_norm()){
    sim_fun = function(n){
        V = covars(n)
        event_time = outcome(covars = V)
        cens_time = cens(covars = V)
        dd = data.table(V)
        names(dd) = paste0("X", 1:ncol(dd))
        dd[, ":=" (true_time = event_time, cens_time = cens_time)]
        dd[, time := pmin(true_time, cens_time)]
        dd[, status := 1*(true_time <= time)]
        return(dd[])
    }
    return(sim_fun)
}
## ----------------- Simulator functions -----------------
## The actual simulators.
## Simulate baseline covariates:
sim_covars_truncated_norm <- function(p = 3, sd = .8, mean = 0, trunc = 1, shift = 0){
    sim_fun = function(n){
        V = matrix(rnorm(n*p, mean = mean, sd = sd), nrow = n, ncol = p)
        V[abs(V)>trunc] = trunc
        V = V+shift
        return(V)
    }
    return(sim_fun)
}
## Functions for simulating conditional event times given covariates
## Simulate data from a model with Weibul as baseline... Using this:
## https://cran.r-project.org/web/packages/survMS/vignettes/how-to-simulate-survival-models.html
sim_coxWeibull <- function(beta = 0, base_shape = 0.7, base_scale = 0.3){
    sim_fun = function(covars){
        nsim = nrow(covars)
        p = ncol(covars)
        if(p>length(beta)) ## Pad with trailing zeros
            beta = c(beta, rep(0, p-length(beta)))
        time = (base_scale*(-log(runif(nsim))/exp(covars %*% beta)))^(1/base_shape)
        return(time)
    }
    return(sim_fun)
}
sim_coxWeibull_quad <- function(beta, base_shape = 0.7, base_scale = 0.3){
    ## beta is now coef of the first component, i.e., beta_1*V1 + beta_2*V1^2 + 0 ...
    sim_fun = function(covars){
        nsim = nrow(covars)
        add_covars = cbind(covars, covars[, 1]^2)
        add_beta = c(beta[1], rep(0, ncol(covars)-1), beta[2])
        time = (base_scale*(-log(runif(nsim))/exp(add_covars %*% add_beta)))^(1/base_shape)
        return(time)
    }
    return(sim_fun)
}
sim_lnorm <- function(beta = 0, sdlog = 1, base_meanlog = 0){
    sim_fun = function(covars){
        nsim = nrow(covars)
        p = ncol(covars)
        if(p>length(beta)) ## Pad with trailing zeros
            beta = c(beta, rep(0, p-length(beta)))
        time = rlnorm(nsim, meanlog = base_meanlog + covars %*% beta, sdlog = sdlog)
        return(time)
    }
    return(sim_fun)
}
sim_lnorm_quad <- function(beta, sdlog = 1){
    ## beta is now coef of the first component, i.e., beta_1*V1 + beta_2*V1^2 + 0 ...
    sim_fun = function(covars){
        nsim = nrow(covars)
        add_covars = cbind(covars, covars[, 1]^2)
        add_beta = c(beta[1], rep(0, ncol(covars)-1), beta[2])
        time = rlnorm(nsim, meanlog = add_covars %*% add_beta, sdlog = sdlog)
        return(time)
    }
    return(sim_fun)
}
## Adaptation from
## https://rdrr.io/github/keaven/nphsim/src/R/rpwexp.r
sim_pwexp <- function(haz_coef, ## dim = ncol(covars) times length(intervals)
                      base_haz = 0.1,
                      intervals=c(0, 3, 6, 9)){
    sim_fun = function(covars){
        all_haz = base_haz*exp(covars %*% haz_coef)
        n_sim = nrow(covars)
        kk = length(intervals)
        tx <- 0
        times <- numeric(n_sim)
        indx <- array(TRUE,n_sim)
        for(ii in 1:kk){
            nindx <- sum(indx)
            if (nindx==0) break
            increment <- rexp(nindx, all_haz[indx, ii])
            times[indx] <- tx + increment
            if (ii<kk){
                tx <- intervals[ii+1]
                indx <- (times > tx)
            }
        }
        return(times)
    }    
    return(sim_fun)
}
## ----------------- visualizer -----------------
visualize_model <- function(simulator, n = 50000, xlim = NULL, n_times = 200, plot = TRUE, show = "both", strat = NULL, ...){
    sim_dat0 = simulator(n = n)
    sim_dat0[, all_true := 1]
    if(is.null(xlim))
        xlim = c(0, sim_dat0[, max(time)])
    time_points = seq(xlim[1], xlim[2], length.out = n_times)
    out_form = Surv(true_time, all_true) ~ 1
    cens_form = Surv(cens_time, all_true) ~ 1
    if(is.null(strat)){
        out_prob = 1-predictRisk(prodlim(out_form, data = sim_dat0), newdata = data.table(1), times = time_points)
        cens_prob = 1-predictRisk(prodlim(cens_form, data = sim_dat0), newdata = data.table(1), times = time_points)
        plot_dat <- rbind(data.table(time = time_points, prob = out_prob, model = "out"),
                          data.table(time = time_points, prob = cens_prob, model = "cens"))
    }else{
        covar_names <- paste(names(strat), collapse = "+")
        out_model <- prodlim(update(out_form, as.formula(paste(".~", covar_names))), data = sim_dat0)
        cens_model <- prodlim(update(cens_form, as.formula(paste(".~", covar_names))), data = sim_dat0)
        plot_dat <- do.call(rbind, lapply(1:nrow(strat), function(strat0){
            out_prob = 1-predictRisk(out_model, newdata = strat[strat0], times = time_points)
            cens_prob = 1-predictRisk(cens_model, newdata = strat[strat0], times = time_points)
            out0 <- cbind(strat[strat0],
                          rbind(data.table(time = time_points, prob = out_prob, model = "out"),
                                data.table(time = time_points, prob = cens_prob, model = "cens")))
            return(out0)
        }))
    }      
    if(!plot){
        return(plot_dat[])
    }else{
        if(is.null(strat))
            col_expr <- expr(NULL)
        else
            col_expr <- expr(factor(get(covar_names)))
        if(show == "both")
            plot_out <- ggplot() + theme_bw() +
                geom_line(data = plot_dat, aes(x = time, y = prob,
                                               col = eval(col_expr)))+
                facet_wrap(~model, ...)
        if(show == "out")
            plot_out <- ggplot() + theme_bw() +
                geom_line(data = plot_dat[model == "out"], aes(x = time, y = prob, col = eval(col_expr)))
        if(show == "cens")
            plot_out <- ggplot() + theme_bw() +
                geom_line(data = plot_dat[model == "cens"], aes(x = time, y = prob, col = eval(col_expr)))
        return(plot_out)
    }    
}


######################################################################
### survival-simulator.R ends here
