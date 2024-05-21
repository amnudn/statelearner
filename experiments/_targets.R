### _targets.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 14 2023 (09:52) 
## Version: 
## Last-Updated: May 21 2024 (09:29) 
##           By: Anders Munch
##     Update #: 45
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
tar_source(here("experiments/zel-sim2-1-n-events-mc-fun.R"))

combine_results <- function(job_name, results_dir){
    fls = list.files(results_dir)
    do.call(rbind,
            lapply(fls[grepl(paste0(job_name, "_taskID_"), fls)],
                   function(fl) fread(paste0("./", results_dir, "/", fl))))
}


list(
    tar_target(results_dir, "results", format = "file"), ## Make it look for changes to results dir
    tar_target(zel_sim2_1, combine_results("zel-sim2-1", results_dir = results_dir)),
    tar_target(zel_sim2_1_n_events, zel_sim2_1_n_events_mc_fun()),
    tar_target(ipcw_fail_sim, {
        run0 = combine_results("ipcw-fail-sim", results_dir = results_dir)
        state_more = combine_results("ipcw-fail-sim-more-state-learners", results_dir = results_dir)
        state_more[, SL := "statelearner_many"]
        rbind(run0, state_more)
    }),
    tar_target(ipcw_fail_sim0_more, combine_results("ipcw-fail-sim-more-state-learners", results_dir = results_dir)),    
    tar_target(ipcw_fail_sim0, combine_results("ipcw-fail-sim", results_dir = results_dir)),
    tar_target(ipcw_fail_sim2, combine_results("ipcw-fail-sim2", results_dir = results_dir)),
    tar_target(ipcw_fail_sim3, combine_results("ipcw-fail-sim3", results_dir = results_dir)),
    tar_target(ipcw_fail_sim_without_rf, combine_results("ipcw-fail-sim-without-rf", results_dir = results_dir))
)


######################################################################
### _targets.R ends here
