### ipcw-fail-sim.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Jul 23 2023 (10:29) 
## Version: 
## Last-Updated: Aug  3 2023 (08:52) 
##           By: Anders Munch
##     Update #: 76
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
library(lava)
library(survival)

sim_data <- function(n){
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

sim_ipcw_fail <- function(ns, eval_time = 20, B = 1, split.method = "cv5"){
    do.call(rbind, lapply(ns, function(nn){
        sim_dt = sim_data(n = nn)
        mms = list(dgm = coxph(Surv(time, event)~strata(A), data = sim_dt, x = TRUE, y = TRUE),
                   km = coxph(Surv(time, event)~1, data = sim_dt, x = TRUE, y = TRUE))
        forms = list(km = Hist(time,event)~1,
                     dgm = Hist(time,event)~strat(A))
        ipcw_briers = lapply(seq_along(forms), function(ii){
            try_message <- try({
                ## What we want to happen
                xx = Score(mms,data = sim_dt,formula = forms[[ii]], metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
                xx = xx[, .(out_model = model, cens_model = names(forms)[ii], loss = Brier, loss_type = "brier", time = eval_time)]
                
            }, silent=TRUE)
            if("try-error"%in%class(try_message)){
                ## What to do if error
                xx = data.table(out_model = names(mms),
                                cens_model = names(forms)[ii],
                                loss = as.numeric(NA),
                                loss_type = "brier",
                                time = eval_time)
            }
            return(xx)
        })                    
        sl = statelearner(learners = list(state = mms,
                                          censoring = list(km = coxph(Surv(time,!event)~1, data = sim_dt,x = TRUE,y = TRUE),
                                                           dgm = coxph(Surv(time,!event)~strata(A), data = sim_dt,x = TRUE,y = TRUE))),
                          data = sim_dt,
                          times = eval_time,
                          B = B,
                          split.method = split.method,
                          integrate = FALSE)
        sl[, ":="(b = NULL, loss_type = "sl-brier", time = eval_time)]
        out = do.call(rbind, c(list(sl), ipcw_briers))        
        out[, n := nn]
        ## Calculate oracle for the two models
        mc_data <- sim_data(n = 20000)
        try_message <- try({
            ## What we want to happen
            oracle = Score(mms, data = mc_data, formula = Hist(eventtime, dummy)~1, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE)$Brier$score
            oracle = oracle[, .(out_model = model, oracle_loss = Brier)]
                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            oracle = data.table(out_model = names(mms), oracle_loss = as.numeric(NA))
        }
        out = merge(out, oracle, by = "out_model", all.x = TRUE)
        return(out[])
    }))
}


######################################################################
### ipcw-fail-sim.R ends here
