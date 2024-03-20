### ipcw-fail-sim.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 15 2023 (16:46) 
## Version: 
## Last-Updated: Jan  5 2024 (15:33) 
##           By: Anders Munch
##     Update #: 29
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
    eval_times = c(10, 20)
    ns = c(300, 600, 1200, 2400)
    ps = c(5)
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
    out = raw_calc[!(type == "cens" & grepl("ipcw", SL))]
    return(out)
}

tt <- sim_fun(1232)


sim_fun <- function(seed){
    ## Write sim setup here
    eval_times = c(18, 36)
    ns = c(600, 1200)
    sim_sets = list(original = simZelefsky_wrapper,
                    indep_cens = simZelefsky_indep_cens_wrapper)
    set.seed(seed)
    raw_calc = do.call(rbind, lapply(seq_along(sim_sets), function(pp){
        sim_set_name = names(sim_sets)[pp]
        dat_sim = function(n) sim_sets[[pp]](n = n)
        true_test = dat_sim(10000)[, !c("time", "status")]
        dummy_true_test = copy(true_test)[, dummy := 1]
        ## full_form0 = formula(paste("Surv(time,status)~", paste(paste0("X", 1:10), collapse = "+")))
        do.call(rbind, lapply(ns, function(nn){
            train0 <- dat_sim(nn)[, !c("true_time", "cens_time")]
            eval0 = eval_sl(train0, true_test, eval_time = eval_times)
            eval0[, ":="(n_obs = nn, sim_set = sim_set_name, seed = seed)]
            ## Get Brier from null model
            out_null = prodlim(Hist(time, status)~1, data = train0)
            cens_null = prodlim(Hist(time, !status)~1, data = train0)
            out_sc = Score(list(out_null), data = dummy_true_test, formula = Hist(true_time, dummy)~1, times = eval_times,null.model = FALSE)$Brier$score
            cens_sc = Score(list(cens_null), data = dummy_true_test, formula = Hist(cens_time, dummy)~1, times = eval_times,null.model = FALSE)$Brier$score
            sc0 = rbind(out_sc[, .(time = times, null_brier = Brier, type = "event")],
                        cens_sc[, .(time = times, null_brier = Brier, type = "cens")])            
            out0 = merge(eval0, sc0, by = c("time", "type"), all.x = TRUE)
            out0[, IPA := 1-brier/null_brier]
            return(out0)
        }))
    }))
    out = raw_calc[!(type == "cens" & grepl("ipcw", SL))] 
    return(out[])
}

tt <- sim_fun(13213)


# Save the results for this task as an individual file in the output folder
fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))


######################################################################
### ipcw-fail-sim.R ends here
