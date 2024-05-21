### _targets.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 14 2023 (09:52) 
## Version: 
## Last-Updated: May 21 2024 (11:38) 
##           By: Anders Munch
##     Update #: 54
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(targets)
library(data.table)
library(here)
tar_source(here("R-code/functions"))

combine_results <- function(job_name, results_dir){
    fls = list.files(results_dir)
    do.call(rbind,
            lapply(fls[grepl(paste0(job_name, "_taskID_"), fls)],
                   function(fl) fread(paste0("./", results_dir, "/", fl))))
}

list(
    tar_target(results_dir, "results", format = "file"), ## Make it look for changes to results dir
    tar_target(zel_sim2_1, combine_results("zel-sim2-1", results_dir = results_dir)),
    tar_target(zel_sim2_1_n_events, {
        ## Get event and censoring probabilities:
        eval_times = seq(6, 36, 6)
        sim_sets = list(original = simZelefsky_wrapper,
                        indep_cens = simZelefsky_indep_cens_wrapper)
        out = do.call(rbind, lapply(seq_along(sim_sets), function(pp){
            sim_set_name = names(sim_sets)[pp]
            mc_sim = sim_sets[[pp]](n = 100000)
            do.call(rbind, lapply(eval_times, function(tt){
                data.table(sim_setting = sim_set_name,
                           time = tt,
                           true_events = mc_sim[, sum(true_time<tt)],
                           true_cens = mc_sim[, sum(cens_time<tt)],
                           at_risk = mc_sim[, sum(time>tt)])[]
            }))
        }))
        out[]
    })
)


######################################################################
### _targets.R ends here
