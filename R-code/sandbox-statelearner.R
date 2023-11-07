### sandbox-statelearner.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  6 2023 (15:48) 
## Version: 
## Last-Updated: Nov  7 2023 (11:45) 
##           By: Anders Munch
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(ggplot2)
ggplot(sl, aes(x = censor, y = loss, shape = cause2)) +
  geom_point(position=position_dodge(width=1), size=.8) +
  theme_bw() + ylab("Integrated Brier score") +
  theme(legend.position="top") +
  xlab("Censoring learner") +
  facet_grid( ~ paste("Cause1 learner =", cause1)) +
  scale_shape_manual("Cause2 learner", values = c(1:4, 6))

sl_B = statelearner(learners = list(cause1 = learners,
                                  cause2 = learners,
                                  censor = learners),
                  data = Melanoma,
                  time = 5*365.25,
                  B = 10)

ggplot(sl_B, aes(x = censor, y = loss, shape = cause2, col = cause2)) +
  geom_errorbar(aes(ymin = loss-2*sd, ymax = loss+2*sd), width = .4,
                position=position_dodge(width=1)) +
  geom_point(position=position_dodge(width=1), size=.8) +
  theme_bw() + ylab("Integrated Brier score") +
  theme(legend.position="top") +
  xlab("Censoring learner") +
  facet_grid( ~ paste("Cause1 learner =", cause1)) +
  scale_shape_manual("Cause2 learner", values = c(1:4, 6))




statelearner(learners = list(cause1 = learners,
                             cause2 = learners,
                             censor = learners),
             data = Melanoma,
             time = 5*365.25,
             B = 5)

strata_cox <- coxph(Surv(time, status == 1)~strata(epicel), data = Melanoma, x = TRUE)
plot(prodlim(Hist(time, 1*(status == 1))~epicel, data = Melanoma))
cox_strata_pred <- predictRisk(strata_cox, newdata = data.table(epicel = c("present", "not present")), times = Melanoma[, sort(time)])
lines(Melanoma[, sort(time)], 1-cox_strata_pred[1,], type = "s", col = "blue", lwd = 3)
lines(Melanoma[, sort(time)], 1-cox_strata_pred[2,], type = "s", col = "red", lwd = 3)

cox <- coxph(Surv(time, status == 1)~epicel, data = Melanoma, x = TRUE)
cox_pred <- predictRisk(cox, newdata = data.table(epicel = c("present", "not present")), times = Melanoma[, sort(time)])
lines(Melanoma[, sort(time)], 1-cox_pred[1,], type = "s", col = "blue", lwd = 3)
lines(Melanoma[, sort(time)], 1-cox_pred[2,], type = "s", col = "red", lwd = 3)
predictCHF(strata_cox, newdata = Melanoma[1:5], times = 1000)



ggplot(ate_est[effect == "ATE"], aes(x = time, y = est)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cause), alpha = 0.3) +
  geom_line(aes(col = cause)) +
  geom_point(aes(col = cause)) +
  geom_hline(yintercept = 0)


wd_dat <- form_data(Melanoma)[, event := NULL]
wd_dat[, A := as.numeric(epicel)-1]
wd_dat[, status := 1*(cause1 == 1) + 2*(cause2 == 1)]

cause1_fit <- rfsrc(update(learners[[sl[1, cause1]]]$x_form, Surv(time, cause1)~.),
                    data = wd_dat,
                    ntree = 50)
## predictCHF(cause1_fit, wd_dat[1:5,2:9], times = c(200,400))
cause2_fit <- coxph(update(learners[[sl[1, cause2]]]$x_form, Surv(time, cause2)~.),
                    data = wd_dat,
                    x = TRUE)
censor_fit <- coxph(update(learners[[sl[1, censor]]]$x_form, Surv(time, censor)~.),
                    data = wd_dat,
                    x = TRUE)
## Fitting treatment mechanism:
treat_fit <- glm(A ~ sex + age, family = binomial, data = wd_dat)

## predictRisk(treat_fit, Melanoma)
## predictTreat(treat_fit, Melanoma)

## Fit with build in function:




######################################################################
### sandbox-statelearner.R ends here
