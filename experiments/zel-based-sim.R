### zel-based-sim.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 14 2023 (09:06) 
## Version: 
## Last-Updated: Nov 15 2023 (13:04) 
##           By: Anders Munch
##     Update #: 36
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(here)
library(data.table)
library(targets)
tar_source(here("R-code/functions"))

sl_sim_est <- function(seed, dat_sim, eval_times, ns = c(300, 600, 1200, 2400)){
    set.seed(seed)
    true_test = dat_sim(10000)[, !c("time", "status")]
    do.call(rbind, lapply(ns, function(nn){
        train0 <- dat_sim(nn)[, !c("true_time", "cens_time")]
        out = eval_sl(train0, true_test, eval_time = eval_times)
        out[, ":="(n_obs = nn, seed = seed)]
        return(out)
    }))
}

zel_sim1 <- function(seed){
    sl_sim_est(seed,
               simZelefsky_wrapper,
               eval_times = seq(1, 36, length.out = 100))
}

######################################################################
### zel-based-sim.R ends here
