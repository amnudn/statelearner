### sandbox-statelearner.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  6 2023 (15:48) 
## Version: 
## Last-Updated: Nov  6 2023 (15:48) 
##           By: Anders Munch
##     Update #: 1
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



######################################################################
### sandbox-statelearner.R ends here
