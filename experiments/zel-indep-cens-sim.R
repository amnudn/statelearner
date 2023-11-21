### zel-indep-cens-sim.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 15 2023 (13:05) 
## Version: 
## Last-Updated: Nov 15 2023 (13:48) 
##           By: Anders Munch
##     Update #: 15
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
    dat_sim = simZelefsky_indep_cens_wrapper
    eval_times = seq(1, 36, length.out = 100)
    ns = c(300, 600, 1200, 2400)
    set.seed(seed)
    true_test = dat_sim(10000)[, !c("time", "status")]
    do.call(rbind, lapply(ns, function(nn){
        train0 <- dat_sim(nn)[, !c("true_time", "cens_time")]
        out = eval_sl(train0, true_test, eval_time = eval_times)
        out[, ":="(n_obs = nn, seed = seed)]
        return(out)
    }))
}

result <- sim_fun(seeds[task_id])

# Save the results for this task as an individual file in the output folder
fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))


######################################################################
### zel-indep-cens-sim.R ends here
