### sandbox.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 14 2023 (11:11) 
## Version: 
## Last-Updated: Nov 20 2023 (17:10) 
##           By: Anders Munch
##     Update #: 86
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(ggplot2)

tar_load(zel_sim1_results, store = here("experiments/_targets"))
tar_load(zel_indep_cens_sim_results, store = here("experiments/_targets"))

time_inc <- zel_sim1_results[, diff(unique(time))[1]]
max_time <- zel_sim1_results[, max(time)[1]]
int_brier_zel_sim1 <- zel_sim1_results[!(type == "cens" & grepl("ipcw", SL)), .(int_brier = sum(100*brier)*time_inc/max_time), .(n_obs, type, SL, seed)]
summ_zel_sim1 <- int_brier_zel_sim1[, .(ave_int_brier = mean(int_brier), se = sd(int_brier)/sqrt(.N)), .(n_obs, type, SL)]
dd_ww <- 0.1

ggplot(summ_zel_sim1,
       aes(x = n_obs, y = ave_int_brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww), aes(ymin = ave_int_brier-1.96*se, ymax = ave_int_brier+1.96*se), width = .05, alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) + geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_wrap(~type, scales = "free_y")

time_inc <- zel_indep_cens_sim_results[, diff(unique(time))[1]]
max_time <- zel_indep_cens_sim_results[, max(time)[1]]
int_brier_zel_indep_cens <- zel_indep_cens_sim_results[!(type == "cens" & grepl("ipcw", SL)), .(int_brier = sum(100*brier)*time_inc/max_time), .(n_obs, type, SL, seed)]
summ_zel_indep_cens <- int_brier_zel_indep_cens[, .(ave_int_brier = mean(int_brier, na.rm = TRUE), se = sd(int_brier, na.rm = TRUE)/sqrt(.N)), .(n_obs, type, SL)]

ggplot(summ_zel_indep_cens,
       aes(x = n_obs, y = ave_int_brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww), aes(ymin = ave_int_brier-1.96*se, ymax = ave_int_brier+1.96*se), width = .05, alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) + geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_wrap(~type, scales = "free_y")

dd_ww <- 0.1
tar_load(ipcw_fail_sim, store = here("experiments/_targets"))
summ_ipcw_fail_sim <- ipcw_fail_sim[, .(ave_scaled_int_brier = mean(scaled_int_brier, na.rm = TRUE), se = sd(scaled_int_brier, na.rm = TRUE)/sqrt(.N)), .(n_obs, n_covar, type, SL)]
ggplot(summ_ipcw_fail_sim[SL != "statelearner_many"],
       aes(x = n_obs, y = ave_scaled_int_brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww),aes(ymin = ave_scaled_int_brier-1.96*se, ymax = ave_scaled_int_brier+1.96*se),width = .1,alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) +
    geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_grid(type~n_covar, scales = "free")

ggplot(summ_ipcw_fail_sim,
       aes(x = n_obs, y = ave_scaled_int_brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww),aes(ymin = ave_scaled_int_brier-1.96*se, ymax = ave_scaled_int_brier+1.96*se),width = .1,alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) +
    geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_grid(type~n_covar, scales = "free")

## Check this out:
tar_load(ipcw_fail_sim_without_rf, store = here("experiments/_targets"))

summ_ipcw_fail_sim_without_rf <- ipcw_fail_sim_without_rf[measure_type == "scaled_int_brier", .(brier = mean(brier, na.rm = TRUE), se = sd(brier, na.rm = TRUE)/sqrt(.N)), .(n_obs, n_covar, type, SL)]
ggplot(summ_ipcw_fail_sim_without_rf,
       aes(x = n_obs, y = brier, col = SL)) +
    theme_bw() +
    geom_errorbar(aes(ymin = brier-1.96*se, ymax = brier+1.96*se), width = .05, alpha = .5) + 
    geom_line() + geom_point() +
    scale_x_continuous(trans='log2') +
    facet_grid(type~n_covar, scales = "free")

tar_load(zel_all_results, store = here("experiments/_targets"))
summ_zel_all_results <- zel_all_results[, .(brier = mean(scaled_int_brier, na.rm = TRUE), se = sd(scaled_int_brier, na.rm = TRUE)/sqrt(.N)), .(n_obs, sim_set, type, SL)]
ggplot(summ_zel_all_results,
       aes(x = n_obs, y = brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww),aes(ymin = brier-1.96*se, ymax = brier+1.96*se), width = .05, alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) + geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_wrap(type~sim_set, scales = "free_y", ncol = 4)

tar_load(zel_all0_results, store = here("experiments/_targets"))
summ_zel_all0_results <- zel_all0_results[, .(brier = mean(scaled_int_brier, na.rm = TRUE), se = sd(scaled_int_brier, na.rm = TRUE)/sqrt(.N)), .(n_obs, sim_set, type, SL)]
ggplot(summ_zel_all0_results,
       aes(x = n_obs, y = brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww),aes(ymin = brier-1.96*se, ymax = brier+1.96*se), width = .05, alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) + geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_wrap(type~sim_set, ncol = 4, scales = "free_y")

tar_load(zel_all2_results, store = here("experiments/_targets"))
summ_zel_all2_results <- zel_all2_results[, .(brier = mean(scaled_int_brier, na.rm = TRUE), se = sd(scaled_int_brier, na.rm = TRUE)/sqrt(.N)), .(n_obs, sim_set, type, SL)]
ggplot(summ_zel_all2_results,
       aes(x = n_obs, y = brier, col = SL)) +
    theme_bw() +
    geom_errorbar(position=position_dodge(width = dd_ww),aes(ymin = brier-1.96*se, ymax = brier+1.96*se), width = .05, alpha = .5) + 
    geom_line(position=position_dodge(width = dd_ww)) + geom_point(position=position_dodge(width = dd_ww)) +
    scale_x_continuous(trans='log2') +
    facet_wrap(type~sim_set, scales = "free_y", ncol = 4)

## Could these things be because we are fitting well to the censoring?
## ... maybe not...
## Could consider variable screening for effect on the outcome model...
## But we also need to fit the censoring wel...
## That is, for targeted learning, it is not clear which of them we fit best...

######################################################################
### sandbox.R ends here
