### cumhazards.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  6 2023 (12:03) 
## Version: 
## Last-Updated: Nov  6 2023 (12:11) 
##           By: Anders Munch
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary:
##
## Function to get cumulative hazard function from fitted objects.
##
## Returns matrix of dimension nrow(newdata) * length(times)
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(prodlim)
library(data.table)

predictCHF <- function(object, newdata, ...){
    UseMethod("predictCHF",object)
}
predictCHF.coxph <- function(object, newdata, times,...){
    pred_obj = predictCox(object, newdata, times = times, type = "cumhazard")
    chf = pred_obj$cumhazard
    ## Fix NA by carrying forward... Maybe not the best...
    last_event_Time_ii = prodlim::sindex(jump.times = pred_obj$times,eval.times = pred_obj$lastEventTime)
    if(last_event_Time_ii<ncol(chf))
        chf[, (last_event_Time_ii+1):ncol(chf)] = chf[, last_event_Time_ii]
    return(chf)
}
predictCHF.GLMnet <- function(object, newdata, times,...){
    ## Hackidy-hack...
    risk_pred = riskRegression:::predictRisk.GLMnet(object = object, newdata = newdata, times = times, product.limit = 0, ...)
    chf = -log(1-risk_pred)
    return(chf)
}
predictCHF.rfsrc <- function(object, newdata, times, ...){
    ## Unsafe hack...
    if("time" %in% names(newdata))
        wd = copy(newdata)[, -c("time"), with = FALSE]
    else
        wd = copy(newdata)
    jump_times = object$time.interest
    cschf = predict(object, newdata = wd)$chf[, , 1]
    if(jump_times[1] != 0){
        jump_times = c(0, jump_times)
        cschf = cbind(0, cschf)
    }
    ii = prodlim::sindex(jump.times = jump_times,eval.times = times)
    return(cschf[, ii])
}

## ## alternative stuff
## construct_pred_fun <- function(object, ...){
##     UseMethod("construct_pred_fun",object)
## }
## construct_pred_fun.coxph <- function(model, ...){
##     out = function(newdata,times)
##         matrix(predictCox(model, newdata = newdata, times = times, type = "cumhazard")$cumhazard, ncol = length(times))
##     return(out)
## }
## construct_pred_fun.glm <- function(model, ...){
##     out = function(newdata)
##         matrix(predict(model, newdata = newdata, type = "response"), ncol = 1)
##     return(out)    
## }


######################################################################
### cumhazards.R ends here
