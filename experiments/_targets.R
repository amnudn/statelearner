### _targets.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 14 2023 (09:52) 
## Version: 
## Last-Updated: Nov 21 2023 (08:33) 
##           By: Anders Munch
##     Update #: 38
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

combine_results <- function(job_name, results_dir){
    fls = list.files(results_dir)
    do.call(rbind,
            lapply(fls[grepl(paste0(job_name, "_taskID_"), fls)],
                   function(fl) fread(paste0("./", results_dir, "/", fl))))
}


list(
    tar_target(results_dir, "results", format = "file"), ## Make it look for changes to results dir
    tar_target(zel_sim1_results, combine_results("zel_sim1", results_dir = results_dir)),
    tar_target(zel_indep_cens_sim_results, combine_results("zel-indep-cens-sim", results_dir = results_dir)),
    tar_target(zel_all_results, combine_results("zel-all-settings", results_dir = results_dir)),
    tar_target(zel_all0_results, combine_results("zel-all-settings0", results_dir = results_dir)),
    tar_target(zel_all2_results, combine_results("zel-all-settings2", results_dir = results_dir)),
    tar_target(ipcw_fail_sim, {
        run0 = combine_results("ipcw-fail-sim", results_dir = results_dir)
        state_more = combine_results("ipcw-fail-sim-more-state-learners", results_dir = results_dir)
        state_more[, SL := "statelearner_many"]
        rbind(run0, state_more)
    }),
    tar_target(ipcw_fail_sim0_more, combine_results("ipcw-fail-sim-more-state-learners", results_dir = results_dir)),    
    tar_target(ipcw_fail_sim0, combine_results("ipcw-fail-sim", results_dir = results_dir)),
    tar_target(ipcw_fail_sim2, combine_results("ipcw-fail-sim2", results_dir = results_dir)),
    tar_target(ipcw_fail_sim_without_rf, combine_results("ipcw-fail-sim-without-rf", results_dir = results_dir))
)


######################################################################
### _targets.R ends here
