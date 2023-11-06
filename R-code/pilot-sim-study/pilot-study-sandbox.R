### pilot-study-sandbox.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  3 2023 (11:50) 
## Version: 
## Last-Updated: Nov  6 2023 (11:59) 
##           By: Anders Munch
##     Update #: 32
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


form1 = Surv(time, event==1)~treat+sex+differ+age+adhere+nodes+nodes.squared+treat*sex+treat*perfor
form2 = Surv(time, event==2)~treat+sex+nodes+differ+age
form0 = Surv(time, event==0)~treat+sex+nodes+treat*sex
formA = treat~sex+age+nodes+differ+nodes.squared
fit.colon.cr <- fit.colon.fun(      
    ## formula.1=Surv(time, event==1)~rx+sex+differ+age+nodes.squared+obstruct+perfor+adhere+extent+surg+rx*sex+rx*perfor+rx*age,
    ## formula.2=Surv(time, event==2)~rx+sex+nodes+differ+age+obstruct+adhere+extent+surg,
    ## formula.0=Surv(time, event==0)~rx+sex+nodes+differ+age+obstruct+perfor+adhere+extent+surg,
    ## formula.treat=rx~sex+age+nodes+differ+obstruct+perfor+adhere+extent+surg,
    formula.1=form1,
    formula.2=form2,
    formula.0=form0,
    formula.treat=formA,
    d=colon_binary)
#-- here we simulate data (fixed seed); 
## set.seed(31)      
sim.colon.cr <- synthesize.colon.fun(fit.colon=fit.colon.cr,
                                     d=colon_binary,
                                     name.treat="treat",
                                     event.name="status",
                                     n = 1000)
## Setup data on right format
sim.colon.cr[, A := as.numeric(treat)-1][, ":="(treat.num = NULL, treat = NULL)]

## Get true stuff:
if(FALSE){
    (true_ate <- synthesize.colon.fun(fit.colon=fit.colon.cr,
                                      get.true.value = 1,
                                      d=colon_binary,
                                      name.treat="treat",
                                      event.name="status",
                                      n = 1000000,
                                      tau = 1000)["F1"]-
         synthesize.colon.fun(fit.colon=fit.colon.cr,
                              get.true.value = 0,
                              d=colon_binary,
                              name.treat="treat",
                              event.name="status",
                              n = 1000000,
                              tau = 1000)["F1"])
}


## fit_form1 <- update(form1, Surv(time, status==1)~factor(A)+.-treat)
## fit_form2 <- update(form2, Surv(time, status==2)~factor(A)+.-treat)
fit_form0 <- update(form0, Surv(time, status==0)~factor(A)+.-treat)
fit_formA <- update(formA, A~.)
## Misspecified
fit_form1 <- update(form1, Surv(time, status==1)~1)
fit_form2 <- update(form2, Surv(time, status==2)~factor(A))
## fit_form0 <- update(form0, Surv(time, status==0)~factor(A)+.-treat)
## fit_form0 <- update(form0, Surv(time, status==0)~1)
## fit_formA <- update(formA, A~1)

## fit_form1 <- update(form1, Surv(time, status==1)~.)
## fit_form2 <- update(form2, Surv(time, status==2)~.)
fit_form0 <- update(form0, Surv(time, status==0)~.)
## fit_formA <- formA

## setup some models:
dd0 <- do.call(colon_simulator, dgm_settings)(1000)

cause1_cox <- coxph(fit_form1,data = dd0,x = TRUE, y = TRUE)
cause2_cox <- coxph(fit_form2,data = dd0,x = TRUE, y = TRUE)
cens_cox <- coxph(fit_form0,data = dd0,x = TRUE, y = TRUE)
treat_model <- glm(fit_formA, family = binomial(), data = dd0)

L1_test <- construct_pred_fun(cause1_cox)
L2_test <- construct_pred_fun(cause2_cox)
L0_test <- construct_pred_fun(cens_cox)
pi_test <- construct_pred_fun(treat_model)

L1_test <- construct_pred_fun(cause1_cox)

L1_test(dd0[1:5], seq(100, 1000, 300))

predictCHF(cause1_cox, dd0[1:5], seq(100, 1000, 300))

## >1 x >1
L1_test(dd0[1:5], times = seq(100, 500, 100))
predictCHF(cause1_cox, dd0[1:5], seq(100, 500, 100))
## > 1 x =1
L1_test(dd0[1:5], times = 100)
predictCHF(cause1_cox, dd0[1:5], 100)
## =1 x >1
L1_test(dd0[1], times = seq(100, 500, 100))
predictCHF(cause1_cox, dd0[1], times = seq(100, 500, 100))
## =1 x =1
L1_test(dd0[1], times = 100)
predictCHF(cause1_cox, dd0[1], times = 100)


## L0_test(sim.colon.cr, times = 1000)
## pi_test(sim.colon.cr)

## sim.colon.cr[status == 0, hist(time)] ## WTF?

raw_est <- one_step(data = sim.colon.cr,
               t = 1000,
               Lambda1 = L1_test,
               Lambda2 = L2_test,
               Gamma = L0_test,
               pi = pi_test,
               jump_points = sort(unique(sim.colon.cr[, time])), collapse = 2)

hh <- abs_risk_ate(data = sim.colon.cr,
                   t = 1000,
                   Lambda1 = L1_test,
                   Lambda2 = L2_test,
                   Gamma = L0_test,
                   pi = pi_test,
                   jump_points = sort(unique(sim.colon.cr[, time])), collapse = 2)

(naiv_est <- mean(raw_est$naiv_i))
(os_est <- mean(raw_est$naiv_i) + mean(raw_est$W_i*raw_est$A_i) - mean(raw_est$W_i*raw_est$B_i) + mean(raw_est$W_i*raw_est$C_i))
## 0.0002660844
## -0.089653

ss <- do.call(colon_simulator, dgm_settings)

dat_test <- ss(1000)

fit_estimator(data = dat_test,
              nuisance_form = nuisance_settings[["mis_treat"]],
              est_times = c(1000, 1500))

tt <- simulate_estimator(n = 1000, m = 1, dgm = dgm_settings, nuisance_form = nuisance_settings,
                         est_times = c(1000,1500))

## mean(raw_est$W_i*raw_est$A_i) - mean(raw_est$W_i*raw_est$B_i) + mean(raw_est$W_i*raw_est$C_i)


## mean(raw_est$A_i)
## mean(raw_est$B_i)
## mean(raw_est$C_i)

## Check bias.
## Check DR property...
## Check that martingale integrals are zero.

## Check estimates
library(ggplot2)
## true_vals[, time := t]
ggplot(sim_est[effect == "ATE"],
       aes(x = n_obs, y = est, fill = type, group = paste(n_obs, type))) +
    geom_boxplot() +
    facet_grid(time~nui_set) +
    theme_bw() +
    geom_hline(data = true_vals[effect == "ATE"], aes(yintercept = true_val), col = "red")

ggplot(sim_est[type == "martingale"],
       aes(x = n_obs, y = est, fill = effect, group = paste(n_obs, effect))) +
    geom_boxplot() +
    facet_grid(time~nui_set) +
    theme_bw() +
    geom_hline(yintercept = 0, col = "red")

coverage_dt <- merge(sim_est[effect == "ATE" & type == "one_step"],
                     true_vals[effect == "ATE", .(time, true_val)],
                     by = "time", all.x = TRUE)[
  , .(coverage = mean(((est-1.96*see) <= true_val) & (true_val <= (est+1.96*see)))),
    .(time, nui_set, n_obs)
]

ggplot(coverage_dt,
       aes(x = n_obs, y = coverage)) +
    geom_line() +
    geom_point() +
    facet_grid(time~nui_set) +
    theme_bw() +
    geom_hline(yintercept = 0.95, col = "red")



######################################################################
### pilot-study-sandbox.R ends here
