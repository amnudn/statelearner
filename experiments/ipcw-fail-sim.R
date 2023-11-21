### ipcw-fail-sim.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 15 2023 (16:46) 
## Version: 
## Last-Updated: Nov 15 2023 (17:10) 
##           By: Anders Munch
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
options(rf.cores = 1) ## Keep RF on one core
setDTthreads(1) ## Avoid DT spawning

## Get input
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
max_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX"))
job_name <- as.character(Sys.getenv("SLURM_JOB_NAME"))

## Get seeds
set.seed(12324)
seeds <- sample(1:1e6, size = max_task_id, replace = FALSE)

## Run sim function with current seed
library(here)
library(targets)
tar_source(here("R-code/functions"))

sim_fun <- function(seed){
    start_t = 1
    end_t = 20
    time_inc = (end_t-start_t)/100
    eval_times = seq(start_t, end_t, time_inc)
    ns = c(300, 600, 1200, 2400)
    ps = c(0, 1, 5, 10)
    set.seed(seed)
    raw_calc = do.call(rbind, lapply(ps, function(pp){
        dat_sim = function(n) ipcw_fail_sim_data(n = n, p = pp)
        true_test = dat_sim(10000)[, !c("time", "status")]
        do.call(rbind, lapply(ns, function(nn){
            train0 <- dat_sim(nn)[, !c("true_time", "cens_time")]
            out = eval_sl(train0, true_test, eval_time = eval_times)
            out[, ":="(n_obs = nn, n_covar = pp, seed = seed)]
            return(out)
        }))
    }))
    out = raw_calc[!(type == "cens" & grepl("ipcw", SL)),
                   .(scaled_int_brier = sum(100*brier)*time_inc/end_t),
                   .(n_obs, n_covar, type, SL, seed)]
    return(out)
}

result <- sim_fun(seeds[task_id])

# Save the results for this task as an individual file in the output folder
fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))


######################################################################
### ipcw-fail-sim.R ends here
