### sandbox-alt-approaches.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  7 2023 (14:30) 
## Version: 
## Last-Updated: Nov 14 2023 (14:11) 
##           By: Anders Munch
##     Update #: 187
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(survSuperLearner)
library(data.table)
library(here)
library(targets)
tar_source(here("R-code/functions/"))

sim_dat <- function(n){
    X <- data.frame(X1 = rnorm(n), X2 = rbinom(n, size = 1, prob = 0.5))     
    S0 <- function(t, x) pexp(t, rate = exp(-2 + x[,1] - x[,2] + .5 * x[,1] * x[,2]), lower.tail = FALSE)
    T <- rexp(n, rate = exp(-2 + X[,1] - X[,2] + .5 *  X[,1] * X[,2]))    
    G0 <- function(t, x) {
        as.numeric(t < 15) * .9 * pexp(t, rate = exp(-2 -.5 * x[,1] - .25 * x[,2] + .5 * x[,1] * x[,2]), lower.tail=FALSE)
    }
    C0 <- rbinom(n, 1, .1)
    C <- rexp(n, exp(-2 -.5 * X[,1] - .25 * X[,2] + .5 * X[,1] * X[,2]))
    ## C[C0 == 1] <- 0 ## This makes no sense?!?
    C[C > 15] <- 15     
    time <- pmin(T, C)
    event <- as.numeric(T <= C)
    setDT(X)
    X[, ":="(time = time, status = event, true_time = T, cens_time = C)]
    return(X[])
}

## Westling approach
dt0 <- sim_dat(200)
event.SL.library <- cens.SL.library <- lapply(c("survSL.km", "survSL.coxph", "survSL.rfsrc"), function(alg) {
    ## c(alg, "survscreen.glmnet", "survscreen.marg", "All")
    c(alg, "All")
})

event.SL.library <- cens.SL.library <- list(
                        km = c("survSL.km", "All"),
                        cox = c("survSL.coxph", "All"),
                        rf = c("survSL.rfsrc", "All"))

## event.SL.library <- cens.SL.library <- lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) {
##     c(alg, "survscreen.glmnet", "survscreen.marg", "All")
## })
dt_new <- sim_dat(1000)

westling_fit <- survSuperLearner(time = dt0[,time],
                                 event = dt0[,status],
                                 X = dt0[, .(X1, X2)],
                                 newX = dt_new[, .(X1, X2)],
                                 new.times = 10, #seq(0, 15, .1),
                                 event.SL.library = event.SL.library,
                                 cens.SL.library = cens.SL.library,
                                 verbose = TRUE)

westling_fit$event.SL.predict

## NB! Convergence issue encountered at some fit...?

names(westling_fit$event.fitLibrary)

## Discrete survSL
names(which.max(westling_fit$event.coef)[1])

predict.survSL.coxph(westling_fit$event.fitLibrary[[names(which.max(westling_fit$event.coef)[1])]],
                     newX = dt_new[1:5],
                     new.times = c(10,20))

predict(westling_fit$event.fitLibrary[[names(which.max(westling_fit$event.coef)[1])]],
        newX = dt_new[1:5],
        new.times = c(10,20))



predictRisk.survSL.coxph <- function(object,newdata,times,...){
    surv = predict.survSL.coxph(object,
                                newX = newdata,
                                new.times = times)
    out = (1-surv)
    rownames(out) = NULL
    return(out)
}
predictRisk.survSL.km <- function(object,newdata,times,...){
    surv = predict.survSL.km(object,
                             newX = newdata,
                             new.times = times)
    out = (1-surv)
    rownames(out) = NULL
    return(out)
}
predictRisk.survSL.rfsrc <- function(object,newdata,times,...){
    surv = predict.survSL.rfsrc(object,
                                newX = newdata,
                                new.times = times)
    out = (1-surv)
    rownames(out) = NULL
    return(out)
}

## Maybe just use the functionality from survSL...

cox0 <- coxph(Surv(time, status)~X1+X2, data = dt0, x = TRUE)

predictRisk(cox0,
            newdata = dt_new[1:5],
            times = c(10,14))

predictRisk(westling_fit$event.fitLibrary[[3]],
            newdata = dt_new[1:5, .(X1, X2)],
            times = c(10,14))

## Statelearner


learners <- list(
    cox = list(model = "cox", x_form = ~X1+X2),
    km = list(model = "cox", x_form = ~1),
    rf = list(model = "rfsrc", x_form = ~X1+X2)
)

learners <- list(
    cox = list(model = "cox"),
    km = list(model = "cox", x_form = ~1),
    rf = list(model = "rfsrc")
)

state_fit = statelearner(learners = list(cause1 = learners,
                                         ## cause2 = learners,
                                         censor = learners),
                         data = dt0,
                         split.method = "cv3",
                         time = 15, B = 1,
                         verbose = TRUE)


train_dt <- sim_dat(200)[, !c("true_time", "cens_time")]
test_dt <- sim_dat(5000)[, !c("time", "status")]

tt <- eval_sl(train_dt, test_dt, eval_time = 1:5)

library(ggplot2)

ggplot(tt, aes(x = time, y = brier, col = SL)) +
    theme_bw() +
    geom_line() + geom_point() +
    facet_wrap(~type)


## data simulator
sim_covars0 <- function(n, pn = 2, pb = 2){
    N = matrix(rnorm(n*pn, mean = 0, sd = 1), nrow = n, ncol = pn)
    B = matrix((rnorm(n*pb, mean = 0, sd = 1) <= 0.35), nrow = n, ncol = pb)
    NB = N*B
    N2 = N^2
    return(cbind(N,B, NB, N2))
}

data_simulator <- build_sim(outcome = sim_coxWeibull(beta = c(0,0,-.8),base_scale = 50),
                            cens = sim_coxWeibull(beta = c(0,0,0.9),base_scale = 50),
                            covars = sim_covars0)

## data_simulator(300)
visualize_model(data_simulator, strat = data.table(X3 = c(0,1)), xlim = c(0,50))

plot(prodlim(Hist(time, event)~A, data = ipcw_sim_data(10000)))
plot(prodlim(Hist(time, event)~A, data = ipcw_sim_data(10000), reverse = TRUE))

## set.seed(3)
train_dt <- data_simulator(200)[, !c("true_time", "cens_time")]
test_dt <- data_simulator(10000)[, !c("time", "status")]
tt <- eval_sl(train_dt, test_dt, eval_time = 10*(1:5))
ggplot(tt, aes(x = time, y = brier, col = SL)) +
    theme_bw() +
    geom_line() + geom_point() +
    facet_wrap(~type)

train_dt <- data_simulator(200)[, !c("true_time", "cens_time")]
test_dt <- data_simulator(5000)[, !c("time", "status")]

system.time(tt <- eval_sl(train_dt, test_dt, eval_time = 10*(1:5)))

ggplot(tt, aes(x = time, y = brier, col = SL)) +
    theme_bw() +
    geom_line() + geom_point() +
    facet_wrap(~type)




library(lava)
ipcw_sim_data <- function(n){
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
    setnames(out, c("eventtime", "censtime", "event", "A"), c("true_time", "cens_time", "status", "X1"))
    return(out[])
}



train_dt <- ipcw_sim_data(500)[, !c("true_time", "cens_time")]
test_dt <- ipcw_sim_data(5000)[, !c("time", "status")]
tt <- eval_sl(train_dt, test_dt, eval_time = seq(5, 20, 5),
              stateL_learners = list(
                  km = list(model = "cox", x_form = ~1),
                  cox = list(model = "cox")
              ),
              SurvSL_learners = list(
                  c("survSL.km", "All"),
                  c("survSL.coxph", "All")
              ),
              ipcw_learners = list(
                  km = coxph(Surv(time,status)~1, data=train_dt, x = TRUE, y = TRUE),
                  cox = coxph(Surv(time,status)~., data=train_dt, x = TRUE, y = TRUE)
              ))



library(ggplot2)

ggplot(sim_zelf_setting[, .(brier = 100*mean(brier), se = 100*sd(brier)/sqrt(.N)), .(time, n_obs, type, SL)],
       aes(x = n_obs, y = brier, col = SL)) +
    theme_bw() +
    ## geom_errorbar(aes(ymin = brier-1.96*se, ymax = brier+1.96*se), width = 10) + 
    geom_line() + geom_point() +
    facet_wrap(type~time, scales = "free_y")

ggplot(sim_zelf_indep_cens[, .(brier = 100*mean(brier,na.rm = TRUE), se = 100*sd(brier, na.rm = TRUE)/sqrt(.N)), .(time, n_obs, type, SL)],
       aes(x = n_obs, y = brier, col = SL)) +
    theme_bw() +
    ## geom_errorbar(aes(ymin = brier-1.96*se, ymax = brier+1.96*se), width = 10) + 
    geom_line() + geom_point() +
    facet_wrap(type~time, scales = "free_y")


tar_source(here("R-code/functions"))
options(rf.cores = 1)

library(pryr)

mem_change({
    train_dt <- simZelefsky_wrapper(2400)[, !c("true_time", "cens_time")]
    test_dt <- simZelefsky_wrapper(10000)[, !c("time", "status")]
    tt <- eval_sl(train_dt, test_dt, eval_time = seq(1, 36, length.out = 100))
})

mem_change(eval_sl(train_dt, test_dt, eval_time = seq(1, 36, length.out = 100)))
## 10 MB

library(peakRAM)
peak_large <- peakRAM({
    train_dt <- simZelefsky_wrapper(2400)[, !c("true_time", "cens_time")]
    test_dt <- simZelefsky_wrapper(5000)[, !c("time", "status")]
    eval_sl(train_dt, test_dt, eval_time = seq(1, 36, length.out = 100))
})

peak_small <- peakRAM({
    train_dt <- simZelefsky_wrapper(600)[, !c("true_time", "cens_time")]
    test_dt <- simZelefsky_wrapper(10000)[, !c("time", "status")]
    eval_sl(train_dt, test_dt, eval_time = seq(1, 36, length.out = 100))})

train_dt <- simZelefsky_wrapper(2400)[, !c("true_time", "cens_time")]
test_dt <- simZelefsky_wrapper(10000)[, !c("time", "status")]
X_test <- test_dt[, !c("true_time", "cens_time")]

peakRAM({
    yy <- statelearner(learners = list(cause1 = list(
                                           km = list(model = "cox", x_form = ~1),
                                           cox = list(model = "cox"),
                                           rfsrc = list(model = "rfsrc.fast", ntree = 50)
                                       ),
                                       censor = list(
                                           km = list(model = "cox", x_form = ~1),
                                           cox = list(model = "cox"),
                                           rfsrc = list(model = "rfsrc.fast", ntree = 50)
                                       )),
                       data = train_dt,
                       split.method = "cv10",
                       time = 36,
                       verbose = FALSE)
    1-predictRisk(yy$fitted_winners$cause1, newdata = X_test, times = seq(1, 36, length.out = 100), cause = 1)
    1-predictRisk(yy$fitted_winners$censor, newdata = X_test, times = seq(1, 36, length.out = 100), cause = 1)
})

peakRAM({
    zz = survSuperLearner(time = train_dt[["time"]],
                          event = train_dt[["status"]],
                          X = train_dt[, !c("time", "status")],
                          newX = X_test,
                          new.times = seq(1, 36, length.out = 5),
                          event.SL.library = list(
                              ## c("survSL.km", "All"),
                              c("survSL.coxph", "All")
                              ## c("survSL.rfsrc_small", "All")
                          ),
                          cens.SL.library = list(
                              ## c("survSL.km", "All"),
                              c("survSL.coxph", "All")
                              ## c("survSL.rfsrc_small", "All")
                          ),
                          verbose = FALSE)
    survSL_event_pred = zz$event.SL.predict
    survSL_cens_pred = zz$cens.SL.predict
})
## Using only coxph as learner and with newX of size n=10.000 uses 1688.9 MiB ram...


source(here("R-code/sandbox-eval-sl-prof.R"))
prof <- lineprof(eval_sl_prof_fun(train_dt, test_dt))
shine(prof)

######################################################################
### sandbox-alt-approaches.R ends here
