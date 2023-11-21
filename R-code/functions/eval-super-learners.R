### eval-super-learners.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 10 2023 (14:28) 
## Version: 
## Last-Updated: Nov 20 2023 (13:54) 
##           By: Anders Munch
##     Update #: 283
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(survSuperLearner)

## Make small rfsrc sl for SurvSL... (must be a better way...)
## survSL.rfsrc_small <- function(time, event, X, newX, new.times, obsWeights, id, ...)
##     survSL.rfsrc(time = time, event = event, X = X, newX = newX, new.times = new.times, obsWeights = obsWeights, id = id, ntree = 50, ...)
survSL.rfsrc_small <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
    data <- data.frame(time, event)
    data <- cbind(data, X)
    fit.rfsrc <- randomForestSRC::rfsrc.fast(Surv(time, event) ~ ., data=data, importance = FALSE, case.wt = obsWeights, forest = TRUE, ntree = 50, ...)
    survs <- predict(fit.rfsrc, newdata=newX, importance='none')$survival
    pred <- t(sapply(1:nrow(survs), function(i) {
        stats::approx(c(0,fit.rfsrc$time.interest), c(1,survs[i,]), method='constant', xout = new.times, rule = 2)$y
    }))

  fit <- list(object = fit.rfsrc)
  class(fit) <- c("survSL.rfsrc")
  out <- list(pred = pred, fit = fit)
  return(out)
}
predict.survSL.rfsrc_small <- function(object, newX, new.times, ...) 
    predict.survSL.rfsrc(object = object, newX = newX, new.times = new.times, ...)

survSL.elastic <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
    data <- data.table(time = time, event = event)
    data <- cbind(data, X)
    form0 <- formula(paste0("Surv(time, event)~", paste(names(X), collapse = "+")))    
    ## fit.GLMnet <- GLMnet(Surv(time, event) ~ ., data=data, alpha = .5, weights = obsWeights, ...)
    fit.GLMnet <- GLMnet(form0, data=data, alpha = .5, weights = obsWeights, ...)
    pred <- 1-predictRisk(fit.GLMnet, newdata = newX, times = new.times)
    fit <- list(object = fit.GLMnet)
    class(fit) <- c("survSL.GLMnet")
    out <- list(pred = pred, fit = fit)
    return(out)
}
predict.survSL.elastic <- function(object, newX, new.times, ...)
    1-predictRisk(object$object, newdata = newX, times = new.times)

survSL.lasso <- function(time, event, X, newX, new.times, obsWeights, id, ...) {
    data <- data.table(time = time, event = event)
    data <- cbind(data, X)
    form0 <- formula(paste0("Surv(time, event)~", paste(names(X), collapse = "+")))    
    ## fit.GLMnet <- GLMnet(Surv(time, event) ~ ., data=data, alpha = .5, weights = obsWeights, ...)
    fit.GLMnet <- GLMnet(form0, data=data, weights = obsWeights, ...)
    pred <- 1-predictRisk(fit.GLMnet, newdata = newX, times = new.times)
    fit <- list(object = fit.GLMnet)
    class(fit) <- c("survSL.GLMnet")
    out <- list(pred = pred, fit = fit)
    return(out)
}
predict.survSL.lasso <- function(object, newX, new.times, ...)
    1-predictRisk(object$object, newdata = newX, times = new.times)

eval_sl <- function(train_data,
                    uncens_test_data,
                    eval_times,
                    stateL_learners = list(
                        km = list(model = "cox", x_form = ~1),
                        cox = list(model = "cox"),
                        rfsrc = list(model = "rfsrc.fast", ntree = 50)
                    ),
                    SurvSL_learners = list(
                        c("survSL.km", "All"),
                        c("survSL.coxph", "All"),
                        c("survSL.rfsrc_small", "All")),
                    ipcw_learners = list(
                        km = coxph(Surv(time,status)~1, data=train_data, x = TRUE, y = TRUE),
                        cox = coxph(Surv(time,status)~., data=train_data, x = TRUE, y = TRUE),
                        rfsrc = rfsrc.fast(Surv(time,status)~., data=train_data, forest = TRUE, ntree = 50)
                    )){
    ## Assumes data is a data.table on the form (time, status, X1, ..., Xd).
    ## wd_train = copy(train_data)
    wd_train = train_data
    X_train = wd_train[, !c("time", "status")]
    time_train = wd_train[, time]
    time_status = wd_train[, status]
    ## Only include covariates, uncensored time and "uncensored" censor time in the test data
    ## X_test = copy(uncens_test_data)[, !c("true_time", "cens_time")]
    X_test = uncens_test_data[, !c("true_time", "cens_time")]
    event_ind = at_risk_fun(uncens_test_data, eval_times, "true_time")
    cens_ind = at_risk_fun(uncens_test_data, eval_times, "cens_time")
    ## Fit SurvSL from Westling et als and eval performance
    inner_eval_fun = function(surv_pred, cens_pred){
        out = data.table(time = rep(eval_times,2),
                         type = rep(c("event", "cens"), each = length(eval_times)),
                         brier = as.numeric(NA))
        if(!is.null(surv_pred)){
            eval_event = colSums((event_ind-surv_pred)^2)/nrow(uncens_test_data)
            out[type == "event", brier := eval_event]
        }
        if(!is.null(cens_pred)){
            eval_cens = colSums((cens_ind-cens_pred)^2)/nrow(uncens_test_data)
            out[type == "cens", brier := eval_cens]
        }
        return(out[])
    }
    if(is.null(SurvSL_learners)){
        eval_survSL = NULL
    }else{
        try_message = try({
            ## What we want to happen
            survSL_fit = survSuperLearner(time = time_train,
                                          event = time_status,
                                          X = X_train,
                                          newX = X_test,
                                          new.times = eval_times,
                                          event.SL.library = SurvSL_learners,
                                          cens.SL.library = SurvSL_learners,
                                          verbose = FALSE)
            survSL_event_pred = survSL_fit$event.SL.predict
            survSL_cens_pred = survSL_fit$cens.SL.predict
            rm(survSL_fit)
            eval_survSL = inner_eval_fun(survSL_event_pred, survSL_cens_pred)
            eval_survSL[, SL := "survSL"]
            ## Add winner (selected / highest ranked model):
            ## [TODO]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            eval_survSL = inner_eval_fun(NULL, NULL)
            eval_survSL[, SL := "survSL"]
        }
    }
    ## Fit state learner
    if(is.null(stateL_learners)){
        eval_stateL = NULL
    }else{
        try_message = try({
            ## What we want to happen
            stateL_fit = statelearner(learners = list(cause1 = stateL_learners,
                                                      censor = stateL_learners),
                                      data = wd_train,
                                      split.method = "cv10",
                                      time = max(eval_times),
                                      verbose = FALSE)
            stateL_event_pred = 1-predictRisk(stateL_fit$fitted_winners$cause1, newdata = X_test, times = eval_times, cause = 1)
            stateL_cens_pred = 1-predictRisk(stateL_fit$fitted_winners$censor, newdata = X_test, times = eval_times, cause = 1)
            rm(stateL_fit)
            eval_stateL = inner_eval_fun(stateL_event_pred, stateL_cens_pred)
            eval_stateL[, SL := "statelearner"]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            eval_stateL = inner_eval_fun(NULL, NULL)
            eval_stateL[, SL := "statelearner"]
        }
    }
    if(is.null(ipcw_learners)){
        eval_ipcw_km = NULL
        eval_ipcw_cox = NULL
        eval_oracle = NULL
    }else{
        ## Fit IPCW with km censoring model
        try_message = try({
            ## What we want to happen
            ipcw_ibs_km = Score(ipcw_learners,
                                data = wd_train,
                                formula = Hist(time,status)~1,
                                metrics = "brier",
                                se.fit = FALSE,
                                times = seq(0, max(eval_times), length.out = 100),
                                contrasts = FALSE,
                                null.model = FALSE,
                                B = 1,
                                summary = "ibs",
                                split.method = "cv10",
                                progress.bar = NULL)$Brier$score[times == max(times)]
                                ipcw_km_event_pred = 1-predictRisk(ipcw_learners[[ipcw_ibs_km[which.min(IBS)][1, model]]], newdata = X_test, times = eval_times)
                                eval_ipcw_km = inner_eval_fun(surv_pred = ipcw_km_event_pred, NULL)
                                eval_ipcw_km[, SL := "ipcw_km"]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            eval_ipcw_km = inner_eval_fun(NULL, NULL)
            eval_ipcw_km[, SL := "ipcw_km"]
        }
        ## Fit IPCW with cox censoring model
        try_message = try({
            ## What we want to happen
            cox_form = formula(paste0("Hist(time,status)~", paste(names(wd_train)[grepl("X", names(wd_train))], collapse = "+")))
            ipcw_ibs_cox = Score(ipcw_learners,
                                 data = wd_train,
                                 formula = cox_form,
                                 metrics = "brier",
                                 se.fit = FALSE,
                                 times = seq(0, max(eval_times), length.out = 100),
                                 contrasts = FALSE,
                                 null.model = FALSE,
                                 B = 1,
                                 summary = "ibs",
                                 split.method = "cv10",
                                 progress.bar = NULL)$Brier$score[times == max(times)]
                                 ipcw_cox_event_pred = 1-predictRisk(ipcw_learners[[ipcw_ibs_cox[which.min(IBS)][1, model]]], newdata = X_test, times = eval_times)
                                 eval_ipcw_cox = inner_eval_fun(surv_pred = ipcw_cox_event_pred, NULL)
                                 eval_ipcw_cox[, SL := "ipcw_cox"]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            eval_ipcw_cox = inner_eval_fun(NULL, NULL)
            eval_ipcw_cox[, SL := "ipcw_cox"]
        }
        ## Calculate oracle
        try_message = try({
            ## What we want to happen
            cens_data = copy(train_data)
            cens_data[, status := 1*!status]
            all_model_fit = do.call(rbind, lapply(seq_along(ipcw_learners), function(ii){
                mm = ipcw_learners[[ii]]
                mm_pred = 1-predictRisk(mm, newdata = X_test, times = eval_times)
                ## Not sure this is OK?
                cc = update(mm, data = cens_data)
                cc_pred = 1-predictRisk(cc, newdata = X_test, times = eval_times)
                all_perf = inner_eval_fun(surv_pred = mm_pred, cc_pred)
                all_perf[, mm_ind := ii]
            }))
            eval_oracle = all_model_fit[, .(brier = min(brier)[1], SL = "oracle"), .(time, type)]
        }, silent=FALSE)
        if("try-error"%in%class(try_message)){
            ## What to do if error
            eval_oracle = inner_eval_fun(NULL, NULL)
            eval_oracle[, SL := "oracle"]
        }
    }
    return(rbind(eval_survSL, eval_stateL, eval_ipcw_km, eval_ipcw_cox, eval_oracle)[])
}



######################################################################
### eval-super-learners.R ends here
