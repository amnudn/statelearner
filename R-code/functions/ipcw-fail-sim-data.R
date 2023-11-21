### ipcw-fail-sim-data.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 15 2023 (16:30) 
## Version: 
## Last-Updated: Nov 15 2023 (16:45) 
##           By: Anders Munch
##     Update #: 22
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
library(lava)

ipcw_fail_sim_data <- function(n, p = 0){
    m <- lvm()
    formula_event <- ~f(X1,-.8)
    formula_cens <- ~f(X1,.9)
    event_scale <- 1/500
    cens_scale <- 1/500
    lava::distribution(m,~X1) <- lava::binomial.lvm(p = .3)
    lava::distribution(m,~censtime) <- lava::coxWeibull.lvm(scale=event_scale)
    lava::distribution(m,~eventtime) <- lava::coxWeibull.lvm(scale=cens_scale)
    m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0), "status")
    lava::regression(m) <- stats::update(formula_event, "eventtime~.")
    lava::regression(m) <- stats::update(formula_cens, "censtime~.")
    out <- setDT(sim(m,n))
    setnames(out, c("censtime", "eventtime"), c("cens_time", "true_time"))
    ## add covariates
    if(p>0){
        X2 = data.table(X2 = rnorm(n, mean = 0.5*out[, X1]))
        if(p>1){
            Xp = as.data.table(matrix(rnorm(n*(p-1)), n, p-1))
            names(Xp) = paste0("X", 1+(p:2))
            X = cbind(Xp, X2)
        }else{
            X = X2
        }
        out = cbind(X, out)
    }        
    return(out[])
}

######################################################################
### ipcw-fail-sim-data.R ends here
