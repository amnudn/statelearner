### sandbox.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 16 2023 (13:40) 
## Version: 
## Last-Updated: Nov 18 2023 (12:43) 
##           By: Anders Munch
##     Update #: 230
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(here)
library(targets)
tar_source(here("R-code/functions"))

library(SuperLearner)
library(concrete)

## Testing data
library(riskRegression)
data(Melanoma, package="riskRegression")
setDT(Melanoma)

## Make list of data sets:
pbc <- as.data.table(survival::pbc)[!is.na(trt), ][, A := trt - 1][, c("time", "status", "A", "age", "sex", "albumin")]
dat_list0 <- list(Melanoma=list(data = copy(Melanoma)[,!"event"][, A := as.numeric(ulcer)-1][, ulcer := NULL][],
                                time = 5*365.25,
                                CR=TRUE),
                  pbc = list(data = copy(pbc),
                             time = 1500,
                             CR = TRUE))

target_with_statelearner <- function(dl){
    stopifnot(c("time", "status", "A") %in% names(dl$data))
    dat_stateL = copy(dl$data)
    covars = names(dat_stateL)[!(names(dat_stateL) %in% c("time", "status", "A"))]
    x_form0 = formula(paste("~", paste(c(covars, "A"), collapse = "+")))
    stateL_learners = list(
        km = list(model = "cox", x_form = ~1),
        km_strat = list(model = "cox", x_form = ~strata(A)),
        ## Get stratified cox to work...
        cox = list(model = "cox", x_form = x_form0),
        rfsrc = list(model = "rfsrc", x_form = x_form0, ntree = 200),
        lasso = list(model = "GLMnet", x_form = x_form0),
        elastic = list(model = "GLMnet", x_form = x_form0, alpha = .5)
    )
    stateL_fit = statelearner(learners = list(cause1 = stateL_learners,
                                              cause2 = stateL_learners,
                                              censor = stateL_learners),
                              data = dat_stateL,
                              time = dl$time)    
    treat_sl = SuperLearner(Y = dat_stateL[, A],
                            X = dat_stateL[, covars, with = FALSE],
                            family = binomial(),
                            SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.xgboost"))
    ate_est = os_abs_risk_ate(data = dat_stateL, 
                              eval_times = dl$time,
                              jump_points = dat_stateL[, sort(unique(time))],
                              fit_1 = stateL_fit$fitted_winner$cause1,
                              fit_2 = stateL_fit$fitted_winner$cause2,
                              fit_cens = stateL_fit$fitted_winner$censor,
                              fit_treat = treat_sl)
    ## Missing implement when no CR
    return(ate_est[effect == "ATE"][])
}


target_with_concrete <- function(dl){
    ConcreteArgs = formatArguments(DataTable = dl$data,
                                   EventTime = "time",
                                   EventType = "status",
                                   Treatment = "A",
                                   Intervention = 0:1,
                                   TargetTime = dl$time,
                                   TargetEvent = NULL,
                                   MaxUpdateIter = 250,
                                   CVArg = list(V = 10),
                                   Verbose = FALSE)
    ConcreteEst = doConcrete(ConcreteArgs)
    ConcreteOut = getOutput(ConcreteEst)
    if(dl$CR){
        out_tmle = data.table(cause = c("cause1", "cause2"),
                              time = dl$time,
                              effect = "ATE",
                              est_type = "tmle",
                              est = c(ConcreteOut[Estimator == "tmle" & Event == 1, .SD[Intervention == "A=1"][["Pt Est"]]-.SD[Intervention == "A=0"][["Pt Est"]]],
                                      ConcreteOut[Estimator == "tmle" & Event == 2, .SD[Intervention == "A=1"][["Pt Est"]]-.SD[Intervention == "A=0"][["Pt Est"]]]),
                              see = c(sd(ConcreteEst[["A=1"]]$IC[Event == 1, IC/sqrt(.N)] - ConcreteEst[["A=0"]]$IC[Event == 1, IC/sqrt(.N)]),
                                      sd(ConcreteEst[["A=1"]]$IC[Event == 2, IC/sqrt(.N)] - ConcreteEst[["A=0"]]$IC[Event == 2, IC/sqrt(.N)])))
        out_gcomp = data.table(cause = c("cause1", "cause2"),
                               time = dl$time,
                               effect = "ATE",
                               est_type = "gcomp",
                               est = c(ConcreteOut[Estimator == "gcomp" & Event == 1, .SD[Intervention == "A=1"][["Pt Est"]]-.SD[Intervention == "A=0"][["Pt Est"]]],
                                       ConcreteOut[Estimator == "gcomp" & Event == 2, .SD[Intervention == "A=1"][["Pt Est"]]-.SD[Intervention == "A=0"][["Pt Est"]]]),
                               see = as.numeric(NA))
    }else{
        out_tmle = data.table(time = dl$time,
                              effect = "ATE",
                              est_type = "tmle",
                              est = ConcreteOut[Estimator == "tmle", .SD[Intervention == "A=1"][["Pt Est"]]-.SD[Intervention == "A=0"][["Pt Est"]]],
                              see = sd(ConcreteEst[["A=1"]]$IC[, IC/sqrt(.N)] - ConcreteEst[["A=0"]]$IC[, IC/sqrt(.N)]))
        out_gcomp = data.table(time = dl$time,
                               effect = "ATE",
                               est_type = "gcomp",
                               est = ConcreteOut[Estimator == "gcomp" & Event == 1, .SD[Intervention == "A=1"][["Pt Est"]]-.SD[Intervention == "A=0"][["Pt Est"]]],
                               see = as.numeric(NA))
    }
    out = rbind(out_gcomp, out_tmle)
    out[, ":="(lower = est-1.96*see, upper = est+1.96*see)]
    return(out[])
}

target_with_riskRegression <- function(dl){
    wd = copy(dl$data)
    covars = names(wd)[!(names(wd) %in% c("time", "status", "A"))]
    wd[ , A := factor(A)]
    x_form0 = formula(paste("~", paste(c(covars, "A"), collapse = "+")))
    ff = as.formula(paste("Hist(time, status)", "~", paste(c(covars, "A"), collapse = "+")))
    m.event =  CSC(ff,data=wd)
    m.event$call$formula = eval(ff)
    m.event =  CSC(update(Hist(time,status)~., x_form0),data=wd)
    ## Need this hack...:
    m.event$call$formula = eval(update(Hist(time,status)~., x_form0))
    m.censor =  coxph(update(Surv(time,status==0)~., x_form0),data=wd, x = TRUE, y = TRUE)
    m.treatment =  glm(formula(paste("A~", paste(c(covars), collapse = "+"))),data=wd,family=binomial(link="logit"))
    ateRobust_cause1 = ate(event = m.event,treatment = m.treatment,censor = m.censor,data = wd,
                           times = dl$time,
                           cause = 1,
                           estimator = c("GFORMULA","AIPTW"))
    ateRobust_cause2 = ate(event = m.event,treatment = m.treatment,censor = m.censor,data = wd,
                           times = dl$time,
                           cause = 2,
                           estimator = c("GFORMULA","AIPTW"))
    out = rbind(ateRobust_cause1$diffRisk[, cause := "cause1"], ateRobust_cause2$diffRisk[, cause := "cause2"])
    out = out[, .(time, cause = cause, effect = "ATE", est_type = estimator, est = estimate, see = se, lower = lower, upper = upper)]
    return(out[])
}

set.seed(212)
os1 <- target_with_statelearner(dat_list0[[1]])
set.seed(212)
tmle1 <- target_with_concrete(dat_list0[[1]])
set.seed(212)
ate_fit1 <- target_with_riskRegression(dat_list0[[1]])
os1[est_type == "one-step"]; tmle1[est_type == "tmle"]; ate_fit1[est_type == "AIPTW"]

data_est <- do.call(rbind, lapply(seq_along(dat_list0), function(ii){
    set0 = dat_list0[[ii]]
    data_name = names(dat_list0)[ii]
    os0 <- target_with_statelearner(set0)
    tmle0 <- target_with_concrete(set0)
    ate_fit0 <- target_with_riskRegression(set0)
    out = rbind(os0[est_type == "one-step"], tmle0[est_type == "tmle"], ate_fit0[est_type == "AIPTW"])
    out[, data_set := data_name]
    return(out[])
}))

## Todo:
##   - set cv pars
##   - put into target
##   - find more data set
##   - consider implementing pure survival and comparing in that case

ggplot(data_est,
       aes(x = est_type, y = est)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) +
    geom_point() +
    facet_grid(data_set~cause) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent)+
    geom_hline(yintercept = 0, linetype = 2)  + coord_flip()


plot(prodlim(Hist(time, status)~A, data = dat_list0[[1]]$data), cause = 1)
abline(v = 1826.25)
csc_fit <- CSC(Hist(time, status)~invasion+ici+epicel+thick+sex+age+logthick+A, data = dat_list0[[1]]$data)
dd1 <- copy(dat_list0[[1]]$data)[, A := 1]
dd0 <- copy(dat_list0[[1]]$data)[, A := 0]
mean(predictRisk(csc_fit, dd1, 1826.25, cause = 1)-predictRisk(csc_fit, dd0, 1826.25, cause = 1))
mean(predictRisk(csc_fit, dd1, 1826.25, cause = 2)-predictRisk(csc_fit, dd0, 1826.25, cause = 2))

set.seed(212)
os2 <- target_with_statelearner(dat_list0[[2]])
tmle2 <- target_with_concrete(dat_list0[[2]])
os2; tmle2

## TODO: Compare with CSC!

## Testing!:
set.seed(212)
dat_list0$test <- list(data = copy(dat_list0[[1]]$data)[, status := 1*(status>0)],
                       time = dat_list0[[1]]$time,
                       CR = FALSE)
tmle_surv <- target_with_concrete(dat_list0[["test"]])
target_with_concrete(dat_list0[[1]])


os1; tmle1
    

data <- data[!is.na(trt), ][, trt := trt - 1]
data <- data

ConcreteArgs = formatArguments(DataTable = dat_list0[[1]]$data,
                               EventTime = "time",
                               EventType = "status",
                               Treatment = "A",
                               Intervention = 0:1,
                               TargetTime = 1500,
                               TargetEvent = 1:2,
                               MaxUpdateIter = 250,
                               CVArg = list(V = 10),
                               Verbose = FALSE)

ConcreteEst = doConcrete(ConcreteArgs)
ConcreteOut = getOutput(ConcreteEst)


Melanoma[,A := as.numeric(ulcer)-1]


## Need super learner for treatment model
sl0 <- SuperLearner(Y = Melanoma[, A], X = Melanoma[, .(invasion,ici,epicel,thick,sex,age,logthick)],
                    SL.library = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.xgboost"))

## Setup general code for running stuff on formatted data
library(concrete)
library(survSuperLearner)

## statelearner:

predictTreat(sl0, Melanoma[, .(invasion,ici,epicel,thick,sex,age,logthick)])


######################################################################
### sandbox.R ends here
