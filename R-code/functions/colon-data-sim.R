### validate-estimator.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  2 2023 (12:07) 
## Version: 
## Last-Updated: Nov  3 2023 (15:17) 
##           By: Anders Munch
##     Update #: 196
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(survival)
library(data.table)
library(here)
library(parallel)
source(here("R-code/functions", "hely-colon-sim.R"))

## ## form1=Surv(time, event==1)~sex+differ+age+adhere
## form1=Surv(time, event==1)~rx+sex+differ+age+adhere
## ## form2=Surv(time, event==2)~sex+nodes+differ+age
## form2=Surv(time, event==2)~rx+sex+nodes+differ+age
## form0=Surv(time, event==0)~rx+sex
## formA=rx~sex+age+nodes+differ

fit_estimator <- function(data, nuisance_form, est_times){
    ## Fit nuisance models
    cause1_fit <- coxph(nuisance_form[["formula_1"]],data = data,x = TRUE, y = TRUE)
    cause2_fit <- coxph(nuisance_form[["formula_2"]],data = data,x = TRUE, y = TRUE)
    cens_fit <- coxph(nuisance_form[["formula_cens"]],data = data,x = TRUE, y = TRUE)
    treat_fit <- glm(nuisance_form[["formula_treat"]], family = binomial(), data = data)
    ## Setup on correct format
    L1 <- construct_pred_fun(cause1_fit)
    L2 <- construct_pred_fun(cause2_fit)
    L0 <- construct_pred_fun(cens_fit)
    pi <- construct_pred_fun(treat_fit)
    ## Fit target parameters
    out = do.call(rbind, lapply(est_times, function(est_tt){
        tar_est = abs_risk_ate(data = data,
                               t = est_tt,
                               Lambda1 = L1,
                               Lambda2 = L2,
                               Gamma = L0,
                               pi = pi,
                               jump_points = sort(unique(data[, time])), collapse = 0)

        effect_names = c("A=1","A=0","ATE")
        naiv1 = mean(tar_est$naiv1_i)
        naiv0 = mean(tar_est$naiv0_i)
        naiv_dt = data.table(type = "naiv",
                             effect = effect_names,
                             est = c(naiv1, naiv0, naiv1-naiv0),
                             see = as.numeric(NA))
        debias_term1 = with(tar_est, mean(W1_i*(A_i - tar_est$B_i + C_i)))
        debias_term0 = with(tar_est, mean(W0_i*(A_i - tar_est$B_i + C_i)))
        see_1 = with(tar_est, sd(naiv1_i + W1_i*(A_i - tar_est$B_i + C_i))/sqrt(nrow(data)))
        see_0 = with(tar_est, sd(naiv0_i + W0_i*(A_i - tar_est$B_i + C_i))/sqrt(nrow(data)))
        see_ate = with(tar_est, sd(naiv1_i + W1_i*(A_i - tar_est$B_i + C_i) - (naiv0_i + W0_i*(A_i - tar_est$B_i + C_i)))/sqrt(nrow(data)))
        one_step_dt = data.table(type = "one_step",
                                 effect = effect_names,
                                 est = c(naiv1+debias_term1, naiv0+debias_term0, naiv1+debias_term1 - (naiv0+debias_term0)),
                                 see = c(see_1, see_0, see_ate))
        martingale_dt = data.table(type = "martingale",
                                   effect = paste0("term_", letters[1:3]),
                                   est = c(mean(tar_est$A_i), mean(tar_est$B_i),mean(tar_est$C_i)),
                                   see = as.numeric(NA))
        out = rbind(naiv_dt, one_step_dt, martingale_dt)
        out[, time := est_tt]
        return(out)
    }))
    return(out[])
}

simulate_estimator <- function(n, m, dgm, nuisance_forms, est_times, n_cores = 1){
    dat_simulator = do.call(colon_simulator, dgm)
    out = do.call(rbind, mclapply(X = 1:m, mc.cores = n_cores, FUN = function(mm){
        dat_mm = dat_simulator(n)
        sim_mm = do.call(rbind, lapply(1:length(nuisance_forms),function(ii){
            args0 = list(data = dat_mm,
                         nuisance_form = nuisance_forms[[ii]],
                         est_times = est_times)
            ## Add try catch!
            try_message = try({
                ## What we want to happen
                fit0 = do.call(fit_estimator, args0)
                fit0[, nui_set := names(nuisance_forms)[ii]]
                
            }, silent=TRUE)
            if("try-error"%in%class(try_message)){
                ## What to do if error
                fit0 = NULL
            }
            return(fit0)
        }))
        sim_mm[, sim_nr := mm]
    }))
    out[, n_obs := n]
}



######################################################################
### validate-estimator.R ends here
