### ipcw-fail-example.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Jul 21 2023 (16:27) 
## Version: 
## Last-Updated: Jul 23 2023 (11:46) 
##           By: Anders Munch
##     Update #: 34
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
library(prodlim)
library(riskRegression)
library(targets)
library(lava)
library(data.table)
library(riskRegression)
library(survival)
library(ranger)

try(setwd("~/research/SuperVision/Anders/survival-loss/statelearner/empirical-study/"))
try(setwd("~/Documents/phd/survival-loss-function/statelearner/empirical-study/"))
tar_source("functions")

## Trying with binary
sim2 <- function(n){
    m <- lvm()
    formula_event <- ~f(A,-.8)
    formula_cens <- ~f(A,.9)
    event_scale <- 1/500
    cens_scale <- 1/500
    lava::distribution(m,~A) <- lava::binomial.lvm(p = .3)
    lava::distribution(m,~censtime) <- lava::coxWeibull.lvm(scale=event_scale)
    lava::distribution(m,~eventtime) <- lava::coxWeibull.lvm(scale=cens_scale)
    m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0), "event")
    lava::regression(m) <- stats::update(formula_event, "eventtime~.")
    lava::regression(m) <- stats::update(formula_cens, "censtime~.")
    out <- setDT(sim(m,n))
    out[,dummy := 1]
    return(out[])
}

eval_time <- 20
dlarge <- sim2(1000)
par(mfrow = c(1,2))
plot(prodlim(Hist(time,event)~A, data = dlarge), xlim = c(0,eval_time), type = "risk")
plot(prodlim(Hist(time,event)~A, data = dlarge, reverse = TRUE), xlim = c(0,eval_time))
plot(prodlim(Hist(time,event)~1, data = dlarge, reverse = TRUE), xlim = c(0,eval_time), add = TRUE, col = "red")
par(mfrow = c(1,1))
mms <- list(dgm = coxph(Surv(time, event)~strata(A), data = dlarge, x = TRUE, y = TRUE),
            km = coxph(Surv(time, event)~1, data = dlarge, x = TRUE, y = TRUE))
Score(mms,data = dlarge,formula = Hist(eventtime,dummy)~1,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = 1, split.method = "cv5")
Score(mms,data = dlarge,formula = Hist(time,event)~1,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = TRUE, B = 1, split.method = "cv5")
## Score(mms,data = dlarge,formula = Hist(time,event)~A,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = TRUE,, B = 1, split.method = "cv5")
statelearner(learners = list(state = mms,
                             censoring = list(km_cens = coxph(Surv(time,!event)~1, data = dlarge,x = TRUE,y = TRUE),
                                              dgm_cens = coxph(Surv(time,!event)~strata(A), data = dlarge,x = TRUE,y = TRUE))),
             data = dlarge,
             times = eval_time,
             B = 1,
             integrate = FALSE)



######################################################################
### ipcw-fail-example.R ends here
