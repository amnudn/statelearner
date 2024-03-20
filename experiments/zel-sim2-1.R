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
    eval_times = seq(6, 36, 6)
    ns = c(300, 600, 1200, 2400)
    sim_sets = list(original = simZelefsky_wrapper,
                    indep_cens = simZelefsky_indep_cens_wrapper)
    set.seed(seed)
    raw_calc = do.call(rbind, lapply(seq_along(sim_sets), function(pp){
        sim_set_name = names(sim_sets)[pp]
        dat_sim = function(n) sim_sets[[pp]](n = n)
        true_test = dat_sim(10000)[, !c("time", "status")]
        dummy_true_test = copy(true_test)[, dummy := 1]
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

result <- sim_fun(seeds[task_id])
## tt <- sim_fun(13213)

fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))
