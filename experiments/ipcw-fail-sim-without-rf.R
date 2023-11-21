### ipcw-fail-sim-without-rf.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 15 2023 (20:45) 
## Version: 
## Last-Updated: Nov 15 2023 (20:58) 
##           By: Anders Munch
##     Update #: 17
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

## Load libraries here:
library(here)
library(targets)
tar_source(here("R-code/functions"))

sim_fun <- function(seed){
    ## Write sim setup here
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
            out = eval_sl(train0, true_test, eval_time = eval_times,
                          stateL_learners = list(
                              km = list(model = "cox", x_form = ~1),
                              cox = list(model = "cox")),
                          SurvSL_learners = list(
                              c("survSL.km", "All"),
                              c("survSL.coxph", "All")),
                          ipcw_learners = list(
                              km = coxph(Surv(time,status)~1, data=train0, x = TRUE, y = TRUE),
                              cox = coxph(Surv(time,status)~., data=train0, x = TRUE, y = TRUE))
                          )
            out[, ":="(n_obs = nn, n_covar = pp, seed = seed)]
            return(out)
        }))
    }))
    int_brier = raw_calc[!(type == "cens" & grepl("ipcw", SL)),
                         .(measure_type = "scaled_int_brier",
                           brier = sum(100*brier)*time_inc/end_t),
                         .(n_obs, n_covar, type, SL, seed)]
    point_brier = raw_calc[time == max(eval_times) & !(type == "cens" & grepl("ipcw", SL)),
                           .(measure_type = paste0("point_brier_time", max(eval_times)),
                             brier = 100*brier),
                           .(n_obs, n_covar, type, SL, seed)]
    out = rbind(int_brier, point_brier)
    return(out)
}

result <- sim_fun(seeds[task_id])

# Save the results for this task as an individual file in the output folder
fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))


######################################################################
### ipcw-fail-sim-without-rf.R ends here
