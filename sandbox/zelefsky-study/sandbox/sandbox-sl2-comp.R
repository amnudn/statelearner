### sandbox-v2.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Aug 18 2023 (16:18) 
## Version: 
## Last-Updated: Aug 21 2023 (14:06) 
##           By: Anders Munch
##     Update #: 81
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## Get CHF for rf

library(randomForestSRC)
library(data.table)
data(wihs, package = "randomForestSRC")
setDT(wihs)
wihs_fix <- copy(wihs)
wihs_fix[, time := time + abs(rnorm(.N, sd = 0.01))]
## dd <- form_data(wihs_fix)

learners1 <- list(cox =  list("cox"),
                  cox_small =  list("cox", x_form = ~ageatfda + cd4nadir),
                  rf = list("rfsrc", ntree = 50)
                  )
learner_all <- list(cause1 = learners1, cause2 = learners1, censor = learners1)

sl0 <- statelearner2(learner_all, data = wihs_fix, time = 5, integrate = FALSE, verbose = TRUE)
sl0

zelefsky
zelefsky[, time := dmos]
zelefsky[status == 0 & vital == "Dead", status := 2]
use_dat <- copy(zelefsky)[, .(time,status,logPSA,stage,ggtot,sDose,hormones)]
zel_learner <- list(
    ## cox =  list("cox"),
    cox_lasso =  list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones), ## Does not work without explicit formula?
    ## cox_ridge =  list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones, alpha = 0),
    cox_elastic =  list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones, alpha = 0.5),
    cox_strata_stage = list("cox", x_form = ~logPSA+strata(stage)+ggtot+sDose+hormones),
    km =  list("cox", x_form = ~ 1),                    
    rf = list("rfsrc", ntree = 200)
)

zel_sl2 <- statelearner2(list(cause1 = zel_learner, cause2 = zel_learner, censor = zel_learner),
                         data = use_dat,
                         time = 36,
                         integrate = TRUE,
                         verbose = TRUE,
                         B = 1)

library(ggplot2)
ggplot(zel_sl2, aes(x = cause1, y = loss, col = cause2)) +
    geom_point() +
    ## geom_errorbar(aes(ymin = loss-2*sd, ymax = loss+2*sd)) +
    theme_bw() +
    facet_grid( ~ censor)

9

## Score(list(km = coxph(Surv(time,status)~1, data=use_dat, x = TRUE, y = TRUE),
##            cox = coxph(Surv(time,status)~logPSA+stage+ggtot+sDose+hormones, data=use_dat, x = TRUE, y = TRUE),
##            cox_small = coxph(Surv(time,status)~logPSA+stage, data=use_dat, x = TRUE, y = TRUE),
##            rf = rfsrc.fast(Surv(time,status)~logPSA+stage+ggtot+sDose+hormones, data=use_dat, forest=TRUE,ntree = 50)),
##       data = use_dat,formula = Hist(time,status)~1, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = 5, split.method = split.method, progress.bar = NULL)

## Try to profile to see how long model fitting versus Brier score calculations take



wihs[, length(unique(time))] ## NB many ties

wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)

wihs[, test := sample(c(0,1), size = .N, replace = TRUE)]
wihs[, test2 := test+1]

rfsrc(Surv(time, test) ~ ., wihs, nsplit = 3, ntree = 100)
rfsrc(Surv(time, test2) ~ ., wihs, nsplit = 3, ntree = 100)
rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
rfsrc(Surv(time, status == 1) ~ ., wihs, nsplit = 3, ntree = 100) ## Why does this not work???

tt <- predict(wihs.obj, newdata = wihs[c(1,500,900), .(ageatfda, idu, black, cd4nadir)])

names(tt)

tt$chf

plot.competing.risk(wihs.obj)

cif <- wihs.obj$cif.oob
Time <- wihs.obj$time.interest
idu <- wihs$idu
cif.haart <- cbind(apply(cif[,,1][idu == 0,], 2, mean),
                   apply(cif[,,1][idu == 1,], 2, mean))
cif.aids  <- cbind(apply(cif[,,2][idu == 0,], 2, mean),
                   apply(cif[,,2][idu == 1,], 2, mean))


dd <- form_data(wihs)

wihs2 <- copy(wihs[, status := NULL])
dd2 <- form_data(wihs2, status = "test")

## Write template for fitter
library(survival)
library(prodlim)


wihs_fix <- copy(wihs)
wihs_fix[, time := time + abs(rnorm(.N, sd = 0.01))]
dd <- form_data(wihs_fix)

uu1 <- fit_cause_model("cox", x_form = ~ageatfda + idu + black + cd4nadir, data = dd, "cause1")
uu2 <- fit_cause_model("rfrsc", x_form = ~ageatfda + idu + black + cd4nadir, data = dd, "cause2")
uuc <- fit_cause_model("cox", x_form = ~ageatfda + idu + black + cd4nadir, data = dd, "censor")

time_seq <- dd[, sort(c(seq(0.1,0.5,length.out = 70), sort(unique(time)[1:100])))]## seq(0.1,0.5,length.out = 70)
ch1 <- predictCHF(uu1, newdata = dd[1:5], times = time_seq)
ch2 <- predictCHF(uu2, newdata = dd[1:5], times = time_seq)
chc <- predictCHF(uuc, newdata = dd[1:5], times = time_seq)
hh <- abs_risk_from_cschf(ch1, ch2, chc)
hh[[1]][, 50] + hh[[2]][, 50] + hh[[3]][, 50]

## Compare to CSC
csc_test <- CSC(Hist(time, status) ~ageatfda + idu + black + cd4nadir, data = wihs_fix)
mm1 <- fit_cause_model("cox", x_form = ~ageatfda + idu + black + cd4nadir, data = dd, "cause1")
mm2 <- fit_cause_model("cox", x_form = ~ageatfda + idu + black + cd4nadir, data = dd, "cause2")
jt0 <- dd[time<7, c(0,sort(unique(time)))]
ch1 <- predictCHF(mm1, newdata = dd[1:5], times = jt0)
ch2 <- predictCHF(mm2, newdata = dd[1:5], times = jt0)
abs_risk <- abs_risk_from_cschf(ch1, ch2)
predictRisk(csc_test, newdata = dd[1:5], times = 2:6, product.limit = FALSE, cause = 1)
abs_risk[[1]][, prodlim::sindex(jump.times = jt0,eval.times = 2:6)]
predictRisk(csc_test, newdata = dd[1:5], times = 2:6, product.limit = FALSE, cause = 2)
abs_risk[[2]][, prodlim::sindex(jump.times = jt0,eval.times = 2:6)]
## Done!

yy1 <- predictRisk(csc_test, newdata = dd, times = 1:10, product.limit = FALSE, cause = 1)
yy2 <- predictRisk(csc_test, newdata = dd, times = 1:10, product.limit = FALSE, cause = 2)

time_seq <- dd[time<5, sort(unique(time))]
ch1 <- predictCHF(uu1, newdata = dd, times = time_seq)
ch2 <- predictCHF(uu2, newdata = dd, times = time_seq)
chc <- predictCHF(uuc, newdata = dd, times = time_seq)
hh <- abs_risk_from_cschf(ch1, ch2, chc)
## hh[[1]][, 1000:1010] + hh[[2]][, 1000:1010] + hh[[3]][, 1000:1010]
summary(as.numeric(hh[[1]] + hh[[2]] + hh[[3]]))
mean(as.numeric(hh[[1]] + hh[[2]] + hh[[3]])>1)
## Not all are below 1... Should be fixed...

length(as.numeric(hh[[1]] + hh[[2]] + hh[[3]]))

time_seq <- dd[, sort(rep(sort(unique(time)[1:100]),2))]
ch1_u <- predictCHF(uu1, newdata = dd[1:5], times = time_seq)
ch2_u <- predictCHF(uu2, newdata = dd[1:5], times = time_seq)
chc_u <- predictCHF(uuc, newdata = dd[1:5], times = time_seq)
hh_u <- abs_risk_from_cschf(ch1_u, ch2_u, chc_u)
hh_u[[1]][, 2*50] + hh_u[[2]][, 2*50] + hh_u[[3]][, 2*50]

summary(as.numeric(hh_u[[1]] + hh_u[[2]] + hh_u[[3]]))


rr <- fit_cause_model("rfsrc", x_form = ~ageatfda + idu + black + cd4nadir, data = dd, "cause1", ntree = 50)
predictCHF(rr, newdata = dd[1:5], times = c(0,5,8,90))


abs_risk_from_cschf(matrix(rnorm(8), nrow = 2, ncol = 4),
             matrix(rnorm(8), nrow = 2, ncol = 4),
             matrix(rnorm(8), nrow = 2, ncol = 4))

matrix(1, nrow = 2, ncol = 2) + matrix(1, nrow = 2, ncol = 2)




######################################################################
### sandbox-v2.R ends here
