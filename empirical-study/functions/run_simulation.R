### run_simulation.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 18 2023 (10:49) 
## Version: 
## Last-Updated: Jul 18 2023 (15:22) 
##           By: Thomas Alexander Gerds
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
run_simulation <- function(n,nlarge = 100000,time_horizon){
    learn_data = draw_data(n = n)
    outcome_models <- list(cox0 = coxph(Surv(time,event = event)~1, data=learn_data, x = TRUE, y = TRUE),
                           cox1 = coxph(Surv(time,event = event)~X1+X2, data=learn_data, x = TRUE, y = TRUE),
                           cox2 = coxph(Surv(time,event = event)~X6+X7, data=learn_data, x = TRUE, y = TRUE),
                           ## cox_dgm = coxph(Surv(time,event = event)~ X1 + X2 + X5 + X6 + X7 + X8 + X11 + X12 + quad.X3 + quad.X4 + quad.X5 + quad.X6,data = learn_data,x = TRUE,y = TRUE),
                           rf = ranger(Surv(time, event = event) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12, data = learn_data,num.tree = 50))
    cens_models <- list(cox0 = coxph(Surv(time,event = event)~1, data=learn_data, x = TRUE, y = TRUE),
                        cox1 = coxph(Surv(time,event = cens_status)~X1+X2, data=learn_data, x = TRUE, y = TRUE),
                        cox2 = coxph(Surv(time,event = cens_status)~X6+X7, data=learn_data, x = TRUE, y = TRUE),
                        ## cox_dgm = coxph(Surv(time,event = cens_status)~ X1 + X2 + X5 + X6 + X7 + X8 + X11 + X12 + quad.X3 + quad.X4 + quad.X5 + quad.X6,data = learn_data,x = TRUE,y = TRUE),
                        rf = ranger(Surv(time, event = cens_status) ~ X1+X2+X6+X7, data = learn_data,num.tree = 50))
    sl = statelearner(learners = list(state = outcome_models,censoring = cens_models),data = learn_data,B = 4,times = time_horizon)
    ipcw.km = Score(outcome_models,data = learn_data,split.method = "cv10",B = 1,formula = Hist(time, event)~1,se.fit = FALSE,null.model=FALSE,
                    metric = "brier",times = time_horizon)$Brier$score
    ipcw.cox = Score(outcome_models,data = learn_data,split.method = "cv10",B = 1,formula = Hist(time, event)~X1+X5+X7,se.fit = FALSE,null.model=FALSE,
                     metric = "brier",times = time_horizon)$Brier$score
    large_data = draw_data(n = nlarge)
    large_data[,dummy := 1]
    selected_models = list(statelearner = outcome_models[[sl[1]$out_model]],
                           IPCW_KM = outcome_models[[setkey(ipcw.km,Brier)[1]$model]],
                           IPCW_Cox = outcome_models[[setkey(ipcw.cox,Brier)[1]$model]])
    x = Score(selected_models,
              formula = Hist(timecause1, event = dummy)~1,
              data = large_data,
              times = time_horizon,
              metrics = "brier",
              null.model = TRUE,
              se.fit = FALSE)
    x
}


######################################################################
### run_simulation.R ends here
