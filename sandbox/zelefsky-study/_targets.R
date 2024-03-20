### _targets.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2023 (11:06) 
## Version: 
## Last-Updated: Mar 13 2024 (09:48) 
##           By: Anders Munch
##     Update #: 177
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(here)
try(setwd(here("zelefsky-study")))
try(setwd("~/research/SuperVision/Anders/survival-loss/statelearner/zelefsky-study/"))
library(targets)
library(tarchetypes)
library(parallel)
tar_option_set(packages = c("glmnet","lava","foreach","data.table","prodlim","survival","riskRegression","ranger"))
tar_source("functions")
list(
    tar_target(test_data,draw_data(484)),
    tar_target(test_fit,{
        outcome_models <- list(cox0 = coxph(Surv(time,event = event)~1, data=test_data, x = TRUE, y = TRUE),
                               cox1 = coxph(Surv(time,event = event)~X1+X2, data=test_data, x = TRUE, y = TRUE),
                               cox2 = coxph(Surv(time,event = event)~X6+X7, data=test_data, x = TRUE, y = TRUE),
                               cox_dgm = coxph(Surv(time,event = event)~ X1 + X2 + X5 + X6 + X7 + X8 + X11 + X12 + quad.X3 + quad.X4 + quad.X5 + quad.X6,data = test_data,x = TRUE,y = TRUE),
                               rf = ranger(Surv(time, event = event) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12, data = test_data,num.tree = 50)
                               ## rf_try = ranger(Surv(time, event = event) ~ X1+X2+X6+X7, mtry = 4, data = test_data),
                               ## rf_deep = ranger(Surv(time, event = event) ~ X1+X2+X6+X7, mtry = 4, min.node.size = 1, data = test_data)
                               )
        cens_models <- list(cox1 = coxph(Surv(time,event = cens_status)~X1+X2, data=test_data, x = TRUE, y = TRUE),
                            cox2 = coxph(Surv(time,event = cens_status)~X6+X7, data=test_data, x = TRUE, y = TRUE),
                            cox_dgm = coxph(Surv(time,event = cens_status)~ X1 + X2 + X5 + X6 + X7 + X8 + X11 + X12 + quad.X3 + quad.X4 + quad.X5 + quad.X6,data = test_data,x = TRUE,y = TRUE),
                            rf = ranger(Surv(time, event = cens_status) ~ X1+X2+X6+X7, data = test_data,num.tree = 50)
                            )
        statelearner(learners = list(state = outcome_models,censoring = cens_models),
                     data = test_data,
                     B = 4,
                     times = 2)
    }),
    tar_target(test_ipcw, {
        outcome_models <- list(cox0 = coxph(Surv(time,event = event)~1, data=test_data, x = TRUE, y = TRUE),
                               cox1 = coxph(Surv(time,event = event)~X1+X2, data=test_data, x = TRUE, y = TRUE),
                               cox2 = coxph(Surv(time,event = event)~X6+X7, data=test_data, x = TRUE, y = TRUE),
                               cox_dgm = coxph(Surv(time,event = event)~ X1 + X2 + X5 + X6 + X7 + X8 + X11 + X12 + quad.X3 + quad.X4 + quad.X5 + quad.X6,data = test_data,x = TRUE,y = TRUE),
                               rf = ranger(Surv(time, event = event) ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12, data = test_data,num.tree = 50))
        x = Score(outcome_models,data = test_data,formula = Hist(time,event)~1,metrics = "brier",se.fit = FALSE,times = 2,contrasts = FALSE,null.model = TRUE)
        x$Brier$score
    }),
    tar_target(zelefsky,{
        zelefsky = setDT(get(load("./data/zelefsky.rda")))
        zelefsky$status=as.numeric(as.character(factor(zelefsky$recur,levels=c("No","Yes"),labels=c("0","1"))))
        zelefsky$logPSA <- log(zelefsky$psa)
        zelefsky$sDose <- as.vector(scale(zelefsky$dose))
        zelefsky$hormonesYes <- zelefsky$hormones=="Yes"
        for (v in levels(zelefsky$stage)){
            zelefsky[,paste("stage",v,sep="")] <- as.numeric(zelefsky$stage==v)
        }
        zelefsky[]
    }),
    tar_target(simulation_input,{
        # survival hazards
        zWeib.surv <- survreg(formula = Surv(dmos, recur != "No") ~ log(psa) + stage + ggtot + sDose + hormones,data = zelefsky)
        # censoring hazards
        zWeib.cens <- survreg(formula = Surv(dmos, recur == "No") ~ log(psa) + stage + ggtot + sDose + hormones, data = zelefsky)
        list(surv = list(coef = coef(zWeib.surv),
                         scale = zWeib.surv$scale),
             cens = list(coef = coef(zWeib.cens),
                         scale = zWeib.cens$scale),
             zelefsky = zelefsky)
    },packages = "survival"),
    tar_target(simulated_zelefsky_data,{
        d = simZelefsky(n = 364,
                        censoring = TRUE,
                        simulation_input = simulation_input)
    },packages = c("lava","data.table")),
    tar_target(zelefsky_fit,{
        pred.form <- Surv(dmos,status)~logPSA+stage+ggtot+sDose+hormones
        cens.form <- Surv(dmos,status == 0)~logPSA+stage+ggtot+sDose+hormones
        mod_surv <- coxph(pred.form,data = simulated_zelefsky_data)
        mod_cens <- coxph(cens.form,data = simulated_zelefsky_data)
        list("survival" = mod_surv,"censoring" = mod_cens)
    }),
    tar_target(zelefsky_statelearner,{
        simulated_zelefsky_data[, time := dmos]
        simulated_zelefsky_data[, event := status]
        ## pred.form <- Surv(time,event)~logPSA+stage+ggtot+sDose+hormones
        ## cens.form <- Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones
        outcome_models <- list(km = coxph(Surv(time,event)~1, data=simulated_zelefsky_data, x = TRUE, y = TRUE),
                               cox_full = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=simulated_zelefsky_data, x = TRUE, y = TRUE),
                               cox_small = coxph(Surv(time,event)~ggtot+sDose+hormones, data=simulated_zelefsky_data, x = TRUE, y = TRUE),
                               rf = ranger(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data = simulated_zelefsky_data, num.tree = 50))
        cens_models <- list(km = coxph(Surv(time,event == 0)~1, data=simulated_zelefsky_data, x = TRUE, y = TRUE),
                            cox_full = coxph(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, data=simulated_zelefsky_data, x = TRUE, y = TRUE),
                            cox_small = coxph(Surv(time,event == 0)~ggtot+sDose+hormones, data=simulated_zelefsky_data, x = TRUE, y = TRUE),
                            rf = ranger(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, data = simulated_zelefsky_data,num.tree = 50))
        statelearner(learners = list(state = outcome_models,censoring = cens_models),
                     data = simulated_zelefsky_data,
                     B = 4,
                     times = 36)
    }),
    tar_target(zelefsky_statelearner_real_data,{
        zelefsky[, time := dmos]
        zelefsky[, event := status]
        ## pred.form <- Surv(time,event)~logPSA+stage+ggtot+sDose+hormones
        ## cens.form <- Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones
        outcome_models <- list(km = coxph(Surv(time,event)~1, data=zelefsky, x = TRUE, y = TRUE),
                               cox_full = coxph(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data=zelefsky, x = TRUE, y = TRUE),
                               cox_small = coxph(Surv(time,event)~ggtot+sDose+hormones, data=zelefsky, x = TRUE, y = TRUE),
                               rf = ranger(Surv(time,event)~logPSA+stage+ggtot+sDose+hormones, data = zelefsky, num.tree = 50))
        cens_models <- list(km = coxph(Surv(time,event == 0)~1, data=zelefsky, x = TRUE, y = TRUE),
                            cox_full = coxph(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, data=zelefsky, x = TRUE, y = TRUE),
                            cox_small = coxph(Surv(time,event == 0)~ggtot+sDose+hormones, data=zelefsky, x = TRUE, y = TRUE),
                            rf = ranger(Surv(time,event == 0)~logPSA+stage+ggtot+sDose+hormones, data = zelefsky,num.tree = 50))
        statelearner(learners = list(state = outcome_models,censoring = cens_models),
                     data = zelefsky,
                     B = 4,
                     times = 36)
    }),
    tar_target(zelefsky_statelearner_real_data_comp,{
        zelefsky[, time := dmos]
        use_dat = copy(zelefsky)[, .(time,status,logPSA,stage,ggtot,sDose,hormones,vital)]
        use_dat[status == 0 & vital == "Dead", status := 2][, vital := NULL]
        zel_learner <- list(
            cox_lasso =  list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones), ## Does not work without explicit formula?
            cox_elastic =  list("GLMnet", x_form = ~logPSA+stage+ggtot+sDose+hormones, alpha = 0.5),
            cox_strata_stage = list("cox", x_form = ~logPSA+strata(stage)+ggtot+sDose+hormones),
            km =  list("cox", x_form = ~ 1),
            rf = list("rfsrc", ntree = 500)
        )
        statelearner2(list(cause1 = zel_learner,
                           cause2 = zel_learner,
                           censor = zel_learner),
                      data = use_dat,
                      time = 36,
                      integrate = TRUE,
                      verbose = TRUE,
                      B = 5)}),
    tar_target(ipcw_fail,
               do.call(rbind, mclapply(1:200, mc.cores = 6, FUN = function(m) {
                   out = sim_ipcw_fail(c(300, 600, 900, 1200, 1500))                   
                   out[, n_sim := m]
                   return(out[])
               }))
               ),
    tar_target(sim_zel_learners,
               do.call(rbind, mclapply(1:200, mc.cores = 4, FUN = function(m) {
                   out = sim_learners_zelefsky(c(200, 500, 1000, 2000),simulation_input = simulation_input)
                   out[, n_sim := m]
                   return(out[])
               }))
               ),
    tar_target(sim_zel_learners2, ## Run on Rao
               do.call(rbind, mclapply(1:200, mc.cores = 20, FUN = function(m) {
                   out = sim_learners_zelefsky2(c(300, 600, 900, 1200, 1500),simulation_input = simulation_input)
                   out[, n_sim := m]
                   return(out[])
               }))
               ),
    tar_target(sim_zel_learners2_dgm_cens, ## Run on Rao
               do.call(rbind, mclapply(1:200, mc.cores = 20, FUN = function(m) {
                   out = sim_learners_zelefsky2_dgm_cens(c(300, 600, 900, 1200, 1500),simulation_input = simulation_input)
                   out[, n_sim := m]
                   return(out[])
               }))
               )
)


######################################################################
### _targets.R ends here
