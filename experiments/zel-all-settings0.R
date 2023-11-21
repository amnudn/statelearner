library(data.table)
options(rf.cores = 1) ## Keep RF on one core
setDTthreads(1) ## Avoid DT spawning

## Get input
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
max_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX"))
job_name <- as.character(Sys.getenv("SLURM_JOB_NAME"))

## Get seeds
set.seed(284219)
seeds <- sample(1:1e6, size = max_task_id, replace = FALSE)

## Load libraries here:
library(here)
library(targets)
tar_source(here("R-code/functions"))

sim_fun <- function(seed){
    ## Write sim setup here
    start_t = 1
    end_t = 36
    time_inc = (end_t-start_t)/100
    eval_times = seq(start_t, end_t, time_inc)
    ns = c(300, 600, 1200, 2400)
    sim_sets = list(original = simZelefsky_wrapper,
                    indep_cens = simZelefsky_indep_cens_wrapper,
                    simple_effect = simZelefsky_simple_effect_wrapper,
                    noise = simZelefsky_noise_wrapper)
    set.seed(seed)
    raw_calc = do.call(rbind, lapply(seq_along(sim_sets), function(pp){
        sim_set_name = names(sim_sets)[pp]
        dat_sim = function(n) sim_sets[[pp]](n = n)
        true_test = dat_sim(10000)[, !c("time", "status")]
        ## full_form0 = formula(paste("Surv(time,status)~", paste(paste0("X", 1:10), collapse = "+")))
        do.call(rbind, lapply(ns, function(nn){
            train0 <- dat_sim(nn)[, !c("true_time", "cens_time")]
            x_form0 = formula(paste("~",paste(names(train0)[grepl("X", names(train0))], collapse = "+")))
            out0 = eval_sl(train0, true_test, eval_time = eval_times)
            out0[, ":="(n_obs = nn, sim_set = sim_set_name, seed = seed)]
            return(out0)
        }))
    }))
    out = raw_calc[!(type == "cens" & grepl("ipcw", SL)),
                   .(scaled_int_brier = sum(100*brier)*time_inc/end_t),
                   .(n_obs, sim_set, type, SL, seed)]
    return(out[])
}

result <- sim_fun(seeds[task_id])

fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))
