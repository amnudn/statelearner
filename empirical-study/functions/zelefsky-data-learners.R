### zelefsky-data-learners.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Jul 23 2023 (14:00) 
## Version: 
## Last-Updated: Aug  2 2023 (15:52) 
##           By: Anders Munch
##     Update #: 135
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
library(randomForestSRC)
source("./functions/simZelefsky.R")
sim_learners_zelefsky <- function(ns,simulation_input, eval_time = 36, B = 1, split.method = "cv5"){
    do.call(rbind, lapply(ns, function(nn){
        sim_dt = simZelefsky(n = nn,censoring = TRUE,simulation_input = simulation_input)
        sim_dt[, ":="(time = dmos, event = status)]
        out_mms = list(km = coxph(Surv(time,event)~1, data=sim_dt, x = TRUE, y = TRUE),
                       cox_full = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_int = coxph(Surv(time,event)~pspline(logPSA) + sDose + stage + ggtot+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_spline = coxph(Surv(time,event)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_small = coxph(Surv(time,event)~logPSA+sDose+stage, data=sim_dt, x = TRUE, y = TRUE)
                       )
        cens_mms = list(km = coxph(Surv(time,event==0)~1, data=sim_dt, x = TRUE, y = TRUE),
                        cox_full = coxph(Surv(time,event==0)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_int = coxph(Surv(time,event==0)~pspline(logPSA) + sDose + stage + ggtot+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_spline = coxph(Surv(time,event==0)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_small = coxph(Surv(time,event==0)~logPSA+sDose+stage, data=sim_dt, x = TRUE, y = TRUE)
                        )
        ## Fit outcome model with km-ipcw
        try_message = try({
            ## What we want to happen
            ipcw_out = Score(out_mms,data = sim_dt,formula = Hist(time,event)~1, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
            ipcw_out = ipcw_out[, .(out_model = model, cens_model = "pre-KM", loss = Brier, learner = "ipcw-km", time = eval_time)]                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            ipcw_out = data.table(out_model = names(out_mms),
                                  cens_model = "pre-KM",
                                  loss = as.numeric(NA),
                                  learner = "ipcw-km",
                                  time = eval_time)
        }
        ## Fit censoring model with km-ipcw
        try_message = try({
            ## What we want to happen
            ipcw_cens = Score(cens_mms,data = sim_dt,formula = Hist(time,event == 0)~1, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
            ipcw_cens = ipcw_cens[, .(out_model = "pre-KM", cens_model = model, loss = Brier, learner = "ipcw-km", time = eval_time)]                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            ipcw_cens = data.table(out_model = "pre-KM",
                                   cens_model = names(cens_mms),
                                   loss = as.numeric(NA),
                                   learner = "ipcw-km",
                                   time = eval_time)
        }
        ## Fit statelearner
        try_message = try({
            ## What we want to happen
            sl = statelearner(learners = list(state = out_mms,
                                              censoring = cens_mms),
                              data = sim_dt,
                              times = eval_time,
                              B = B,
                              split.method = split.method,
                              integrate = FALSE)
            sl[, ":="(b = NULL, learner = "sl", time = eval_time)]
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            sl = expand.grid(out_model = names(out_mms),
                             cens_model = names(cens_mms))
            setDT(sl)            
            sl[ , ":="(loss = as.numeric(NA), learner = "sl", time = eval_time)]
        }

        ## Calculate oracle
        mc_data <- simZelefsky(n = 20000,censoring = TRUE,simulation_input = simulation_input)
        mc_data[, ":="(time = dmos, event = status)]
        mc_data[, dummy := 1]
        try_message = try({
            ## What we want to happen
            oracle_out = Score(out_mms,data = mc_data,formula = Hist(time.recur,dummy)~1,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE)$Brier$score
            oracle_out = oracle_out[, .(out_model = model, cens_model = "no-cens", loss = Brier, learner = "oracle", time = eval_time)]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            oracle_out = data.table(out_model = names(out_mms),cens_model = "no-cens",loss = as.numeric(NA),learner = "oracle",time = eval_time)
        }
        try_message = try({
            ## What we want to happen
            oracle_cens = Score(cens_mms,data = mc_data,formula = Hist(time.cens,dummy)~1,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE)$Brier$score
            oracle_cens = oracle_cens[, .(out_model = "no-cens",cens_model = model,loss = Brier,learner = "oracle",time = eval_time)]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            oracle_cens = data.table(out_model = "no-cens",cens_model = names(cens_mms),loss = as.numeric(NA),learner = "oracle",time = eval_time)            
        }
        out = do.call(rbind, list(sl, ipcw_out, ipcw_cens, oracle_out, oracle_cens))
        out[, n := nn]
        return(out[])
    }))
}

sim_learners_zelefsky2 <- function(ns,simulation_input, eval_time = 36, B = 1, split.method = "cv5"){
    options(rf.cores = 1)
    do.call(rbind, lapply(ns, function(nn){
        sim_dt = simZelefsky(n = nn,censoring = TRUE,simulation_input = simulation_input)
        sim_dt[, ":="(time = dmos, event = status)]
        out_mms = list(km = coxph(Surv(time,event)~1, data=sim_dt, x = TRUE, y = TRUE),
                       cox_full = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_strata_stage = coxph(Surv(time,event)~logPSA+strata(stage)+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_strata_hormones = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+strata(hormones), data=sim_dt, x = TRUE, y = TRUE),
                       cox_spline = coxph(Surv(time,event)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       rf = rfsrc.fast(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, forest=TRUE,ntree = 50),
                       lasso = GLMnet(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,data=sim_dt),
                       ridge = GLMnet(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,alpha = 0,data=sim_dt),
                       elastic = GLMnet(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5,data=sim_dt))
        cens_mms = list(km = coxph(Surv(time,event==0)~1, data=sim_dt, x = TRUE, y = TRUE),
                        cox_full = coxph(Surv(time,event==0)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_strata_stage = coxph(Surv(time,event == 0)~logPSA+strata(stage)+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_strata_hormones = coxph(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+strata(hormones), data=sim_dt, x = TRUE, y = TRUE),
                        cox_spline = coxph(Surv(time,event==0)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        rf = rfsrc.fast(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, forest = TRUE,ntree = 50),
                        lasso = GLMnet(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,data=sim_dt),
                        ridge = GLMnet(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,alpha = 0,data=sim_dt),
                        elastic = GLMnet(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5,data=sim_dt)
                        )
        ## Fit outcome model with km-ipcw
        try_message = try({
            ## What we want to happen
            ipcw_out = Score(out_mms,data = sim_dt,formula = Hist(time,event)~1, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
            ipcw_out = ipcw_out[, .(out_model = model, cens_model = "pre-KM", loss = Brier, learner = "ipcw-km", time = eval_time)]                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            ipcw_out = data.table(out_model = names(out_mms),
                                  cens_model = "pre-KM",
                                  loss = as.numeric(NA),
                                  learner = "ipcw-km",
                                  time = eval_time)
        }
        ## Fit censoring model with km-ipcw
        try_message = try({
            ## What we want to happen
            ipcw_cens = Score(cens_mms,data = sim_dt,formula = Hist(time,event == 0)~1, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
            ipcw_cens = ipcw_cens[, .(out_model = "pre-KM", cens_model = model, loss = Brier, learner = "ipcw-km", time = eval_time)]                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            ipcw_cens = data.table(out_model = "pre-KM",
                                   cens_model = names(cens_mms),
                                   loss = as.numeric(NA),
                                   learner = "ipcw-km",
                                   time = eval_time)
        }
        ## Fit statelearner
        try_message = try({
            ## What we want to happen
            sl = statelearner(learners = list(state = out_mms,
                                              censoring = cens_mms),
                              data = sim_dt,
                              times = eval_time,
                              B = B,
                              split.method = split.method,
                              integrate = FALSE)
            sl[, ":="(b = NULL, learner = "sl", time = eval_time)]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            sl = expand.grid(out_model = names(out_mms),
                             cens_model = names(cens_mms))
            setDT(sl)            
            sl[ , ":="(loss = as.numeric(NA), learner = "sl", time = eval_time)]
        }

        ## Calculate oracle
        mc_data <- simZelefsky(n = 20000,censoring = TRUE,simulation_input = simulation_input)
        mc_data[, ":="(time = dmos, event = status)]
        mc_data[, dummy := 1]
        try_message = try({
            ## What we want to happen
            oracle_out = Score(out_mms,data = mc_data,formula = Hist(time.recur,dummy)~1,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE)$Brier$score
            oracle_out = oracle_out[, .(out_model = model, cens_model = "no-cens", loss = Brier, learner = "oracle", time = eval_time)]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            oracle_out = data.table(out_model = names(out_mms),cens_model = "no-cens",loss = as.numeric(NA),learner = "oracle",time = eval_time)
        }
        try_message = try({
            ## What we want to happen
            oracle_cens = Score(cens_mms,data = mc_data,formula = Hist(time.cens,dummy)~1,metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE)$Brier$score
            oracle_cens = oracle_cens[, .(out_model = "no-cens",cens_model = model,loss = Brier,learner = "oracle",time = eval_time)]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            oracle_cens = data.table(out_model = "no-cens",cens_model = names(cens_mms),loss = as.numeric(NA),learner = "oracle",time = eval_time)            
        }
        out = do.call(rbind, list(sl, ipcw_out, ipcw_cens, oracle_out, oracle_cens))
        out[, n := nn]
        return(out[])
    }))
}

sim_learners_zelefsky2_dgm_cens <- function(ns,simulation_input, eval_time = 36, B = 1, split.method = "cv5"){
    options(rf.cores = 1)
    do.call(rbind, lapply(ns, function(nn){
        sim_dt = simZelefsky(n = nn,censoring = TRUE,simulation_input = simulation_input)
        sim_dt[, ":="(time = dmos, event = status)]
        out_mms = list(km = coxph(Surv(time,event)~1, data=sim_dt, x = TRUE, y = TRUE),
                       cox_full = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_strata_stage = coxph(Surv(time,event)~logPSA+strata(stage)+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       cox_strata_hormones = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+strata(hormones), data=sim_dt, x = TRUE, y = TRUE),
                       cox_spline = coxph(Surv(time,event)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_dt, x = TRUE, y = TRUE),
                       rf = rfsrc.fast(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, forest=TRUE,ntree = 50),
                       lasso = GLMnet(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,data=sim_dt),
                       ridge = GLMnet(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,alpha = 0,data=sim_dt),
                       elastic = GLMnet(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5,data=sim_dt))
        cens_mms = list(km = coxph(Surv(time,event==0)~1, data=sim_dt, x = TRUE, y = TRUE),
                        cox_full = coxph(Surv(time,event==0)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_strata_stage = coxph(Surv(time,event == 0)~logPSA+strata(stage)+ggtot+sDose+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        cox_strata_hormones = coxph(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+strata(hormones), data=sim_dt, x = TRUE, y = TRUE),
                        cox_spline = coxph(Surv(time,event==0)~pspline(logPSA)+stage+ggtot+pspline(sDose)+hormones, data=sim_dt, x = TRUE, y = TRUE),
                        rf = rfsrc.fast(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, data=sim_dt, forest = TRUE,ntree = 50),
                        lasso = GLMnet(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,data=sim_dt),
                        ridge = GLMnet(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,alpha = 0,data=sim_dt),
                        elastic = GLMnet(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones,alpha = 0.5,data=sim_dt)
                        )
        ## Fit outcome model with cox-ipcw
        try_message = try({
            ## What we want to happen
            ipcw_out = Score(out_mms,data = sim_dt,formula = Hist(time,event)~logPSA+stage+ggtot+sDose+hormones, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
            ipcw_out = ipcw_out[, .(out_model = model, cens_model = "pre-cox", loss = Brier, learner = "ipcw-cox", time = eval_time)]                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            ipcw_out = data.table(out_model = names(out_mms),
                                  cens_model = "pre-cox",
                                  loss = as.numeric(NA),
                                  learner = "ipcw-cox",
                                  time = eval_time)
        }
        ## Fit censoring model with cox-ipcw
        try_message = try({
            ## What we want to happen
            ipcw_cens = Score(cens_mms,data = sim_dt,formula = Hist(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, metrics = "brier",se.fit = FALSE,times = eval_time,contrasts = FALSE,null.model = FALSE, B = B, split.method = split.method, progress.bar = NULL)$Brier$score
            ipcw_cens = ipcw_cens[, .(out_model = "pre-cox", cens_model = model, loss = Brier, learner = "ipcw-cox", time = eval_time)]                
        }, silent=TRUE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            ipcw_cens = data.table(out_model = "pre-cox",
                                   cens_model = names(cens_mms),
                                   loss = as.numeric(NA),
                                   learner = "ipcw-km",
                                   time = eval_time)
        }
        out = do.call(rbind, list(ipcw_out, ipcw_cens))
        out[, n := nn]
        return(out[])
    }))
}

######################################################################
### zelefsky-data-learners.R ends here
