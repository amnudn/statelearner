### _targets.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Mar 13 2024 (09:20) 
## Version: 
## Last-Updated: Mar 20 2024 (08:58) 
##           By: Anders Munch
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(here)
try(setwd(here("zelefsky-case-study")))
try(setwd("~/research/SuperVision/Anders/survival-loss/statelearner/zelefsky-case-study/"))
library(targets)
library(tarchetypes)
library(parallel)
tar_option_set(packages = c("glmnet","lava","foreach","data.table","prodlim","survival","riskRegression","ranger"))
tar_source(here("R-code/functions"))

list(
    tar_target(zelefsky,{
        zelefsky = setDT(get(load(here("data/zelefsky.rda"))))
        zelefsky$status=as.numeric(as.character(factor(zelefsky$recur,levels=c("No","Yes"),labels=c("0","1"))))
        zelefsky$logPSA <- log(zelefsky$psa)
        zelefsky$sDose <- as.vector(scale(zelefsky$dose))
        zelefsky$hormonesYes <- zelefsky$hormones=="Yes"
        for (v in levels(zelefsky$stage)){
            zelefsky[,paste("stage",v,sep="")] <- as.numeric(zelefsky$stage==v)
        }
        zelefsky[]
    }),
    tar_target(use_data,{
        zelefsky[, time := dmos]
        out = copy(zelefsky)[, .(time,status,logPSA,stage,ggtot,sDose,hormones,vital)]
        out[status == 0 & vital == "Dead", status := 2][, vital := NULL]
        out
    }),
    tar_target(zelefsky_statelearner,{       
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
        statelearner(list(cause1 = zel_learner,
                          cause2 = zel_learner,
                          censor = zel_learner),
                     data = use_dat,
                     time = 36,
                     integrate = TRUE,
                     verbose = TRUE,
                     B = 5)}),
    ## Targeted parameter
    tar_target(ate_est_main_eff,{
        use_data[,A := as.numeric(hormones)-1]
        treat_fit = GLMnet(A ~ logPSA+stage+ggtot+sDose, family = binomial, data = use_data)
        os_abs_risk_ate(data = use_data, 
                        eval_times = seq(6, 36, 6),
                        fit_1 = zelefsky_statelearner$fitted_winners$cause1,
                        fit_2 = zelefsky_statelearner$fitted_winners$cause2,
                        fit_cens = zelefsky_statelearner$fitted_winners$censor,
                        fit_treat = treat_fit)}
        ),
    tar_target(ate_est_inter_eff,{
        use_data[,A := as.numeric(hormones)-1]
        treat_fit = GLMnet(A ~ logPSA*stage*ggtot*sDose, family = binomial, data = use_data)
        os_abs_risk_ate(data = use_data, 
                        eval_times = seq(6, 36, 6),
                        fit_1 = zelefsky_statelearner$fitted_winners$cause1,
                        fit_2 = zelefsky_statelearner$fitted_winners$cause2,
                        fit_cens = zelefsky_statelearner$fitted_winners$censor,
                        fit_treat = treat_fit)}
        )
)


######################################################################
### _targets.R ends here
