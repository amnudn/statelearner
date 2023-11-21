### zel-all-settings.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 20 2023 (09:39) 
## Version: 
## Last-Updated: Nov 20 2023 (21:16) 
##           By: Anders Munch
##     Update #: 75
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
            ## Unnice hack...
            if(sim_set_name == "noise")
                ipcw_glm0 = GLMnet(Surv(time, status)~X6+X7+X8+X9+X10+X1+X2+X3+X4+X5, data = train0, alpha = 0.5)
            else
                ipcw_glm0 = GLMnet(Surv(time, status)~X1+X2+X3+X4+X5, data = train0, alpha = 0.5)            
            out0 = eval_sl(train0, true_test, eval_time = eval_times,
                           stateL_learners = list(
                               km = list(model = "cox", x_form = ~1),
                               cox = list(model = "cox", x_form = x_form0),
                               elastic = list(model = "GLMnet", x_form = x_form0, alpha = .5),
                               rfsrc = list(model = "rfsrc.fast", x_form = x_form0, ntree = 50)),
                           SurvSL_learners = list(
                               c("survSL.km", "All"),
                               c("survSL.coxph", "All"),
                               c("survSL.elastic", "All"),
                               c("survSL.rfsrc_small", "All")
                           ),
                           ## stateL_learners = NULL,
                           ## SurvSL_learners = NULL,
                           ipcw_learners = list(
                               km = coxph(Surv(time,status)~1, data=train0, x = TRUE, y = TRUE),
                               cox = coxph(Surv(time,status)~., data=train0, x = TRUE, y = TRUE),
                               elastic = ipcw_glm0,
                               rfsrc = rfsrc.fast(Surv(time,status)~., data=train0, forest = TRUE, ntree = 50)
                           )
                           )
            out0[, ":="(n_obs = nn, sim_set = sim_set_name, seed = seed)]
            return(out0)
        }))
    }))
    out = raw_calc[!(type == "cens" & grepl("ipcw", SL)),
                   .(scaled_int_brier = sum(100*brier)*time_inc/end_t),
                   .(n_obs, sim_set, type, SL, seed)]
    return(out[])
}

## ## Testing:
## library(peakRAM)
## system.time({
##     peak0<-peakRAM({tt<-sim_fun(131)})
## })

result <- sim_fun(seeds[task_id])

# Save the results for this task as an individual file in the output folder
fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))


######################################################################
### zel-all-settings.R ends here
