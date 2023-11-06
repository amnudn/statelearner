### _targets.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  3 2023 (09:35) 
## Version: 
## Last-Updated: Nov  6 2023 (09:06) 
##           By: Anders Munch
##     Update #: 103
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
try(setwd("~/Documents/phd/statelearner/R-code/pilot-sim-study/"))
library(here)
library(targets)
library(tarchetypes)
library(parallel)
tar_source(here("R-code/functions"))

## DGM
dgm_settings <- list(
    formula_1 = Surv(time, event==1)~treat+sex+differ+age+adhere+nodes+nodes.squared+treat*sex+treat*perfor,
    formula_2 = Surv(time, event==2)~treat+sex+nodes+differ+age,
    formula_cens = Surv(time, event==0)~treat+sex+nodes+treat*sex,
    formula_treat = treat~sex+age+nodes+differ+nodes.squared)

## Sim settings
true_settings <- list(
    formula_1 = Surv(time, status==1)~factor(A)+sex+differ+age+adhere+nodes+nodes.squared+factor(A)*sex+factor(A)*perfor,
    formula_2 = Surv(time, status==2)~factor(A)+sex+nodes+differ+age,
    formula_cens = Surv(time, status==0)~factor(A)+sex+nodes+factor(A)*sex,
    formula_treat = A~sex+age+nodes+differ+nodes.squared)
nuisance_settings <- rep(list(true_settings), 5)
## names(nuisance_settings) <- c("all_correct", "mis_treat", "mis_out")
## names(nuisance_settings) <- c("all_correct", "mis_treat", "mis_cens", "mis_out")
names(nuisance_settings) <- c("all_correct", "mis_treat", "mis_cens", "mis_out", "mis_all")
nuisance_settings[["mis_treat"]][["formula_treat"]] <- A~1
nuisance_settings[["mis_cens"]][["formula_cens"]] <- Surv(time, status==0)~1
nuisance_settings[["mis_out"]][["formula_1"]] <- Surv(time, status==1)~factor(A)
nuisance_settings[["mis_out"]][["formula_2"]] <- Surv(time, status==1)~factor(A) + perfor
nuisance_settings[["mis_all"]][["formula_treat"]] <- A~1
nuisance_settings[["mis_all"]][["formula_cens"]] <- Surv(time, status==0)~1
nuisance_settings[["mis_all"]][["formula_1"]] <- Surv(time, status==1)~factor(A)
nuisance_settings[["mis_all"]][["formula_2"]] <- Surv(time, status==1)~factor(A) + perfor

eval_times <- c(1000, 2000)
n_vals <- c(1000, 1500, 2000)

list(
    tar_target(true_vals,{
        dat_sim = do.call(colon_simulator, dgm_settings)
        do.call(rbind,lapply(eval_times, function(tt){
            data.table(effect = c("A=1","A=0","ATE"),
                       true_val = unlist(dat_sim(1000000, get_true_ate = TRUE, t = tt)),
                       time = tt)
        }))
    }),
    ## Set up with tar_rep instead later
    tar_rep(sim_est,
    {
        do.call(rbind, lapply(n_vals, function(nn){
            simulate_estimator(n = nn, m = 24, dgm = dgm_settings, nuisance_forms = nuisance_settings, est_times = eval_times, n_cores = 6)
        }))
    },
    batches = 10)
)


######################################################################
### _targets.R ends here
