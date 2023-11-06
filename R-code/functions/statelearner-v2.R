### statelearner-v2.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Aug 18 2023 (10:26) 
## Version: 
## Last-Updated: Nov  6 2023 (13:03) 
##           By: Anders Munch
##     Update #: 366
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


## Set data on canonical form
form_data <- function(data, max_time = Inf, time = "time", status = "status", cause_codes = c("1" = 1, "2" = 2, "c" = 0), vars = NULL){
    if(!is.null(vars))
        wd = copy(data)[, c(time, status, vars), with = FALSE]
    else
        wd = copy(data)
    setnames(wd, old = c(time, status), new = c("time", "status"))
    ## Make cause of interest = 1, everything else = 2
    wd[, cause1 := 1+1*(status != cause_codes[["1"]] & time <= max_time)]
    wd[, censor := 1+1*(status != cause_codes[["c"]] & time <= max_time)]
    if(wd[, length(unique(status))>2])
        wd[, cause2 := 1+1*(status != cause_codes[["2"]] & time <= max_time)]
    wd[, status := NULL]
    return(wd[])
}
## Fit model on canonical data form
fit_cause_model <- function(model, data, cause, x_form = NULL, ...){
    cause_names = c("cause1", "cause2", "censor")
    cause_names = cause_names[cause_names %in% names(data)]
    stopifnot(cause %in% cause_names)
    wd = copy(data)[, -cause_names[cause_names != cause], with = FALSE]
    setnames(wd, old = cause, new = "event")
    ## not so elegant
    if(is.null(x_form))
        form = Surv(time, event) ~ .
    else
        form = update(x_form, Surv(time, event) ~ .)
    if(model == "cox"){
        wd[, event := 1*(event == 1)] ## To make this cause of interest
        out = coxph(form, data = wd,x = TRUE,y = TRUE, ...)
    }
    if(model == "GLMnet"){
        wd[, event := 1*(event == 1)] ## To make this cause of interest
        out = GLMnet(form, data = wd, ...)
    }
    if(model == "rfsrc"){
        out = rfsrc(form, data = wd, ...)
    }
    return(out)
}
## Calculate F from CSCHF's
abs_risk_from_cschf <- function(...){
    chfs = list(...)
    S = exp(-Reduce("+", chfs))
    Sminus = cbind(1,S[,-ncol(S)])
    abs_risk = lapply(chfs, function(xx) {
        cs_haz = t(apply(cbind(0,xx), 1, diff))
        abs_risk_diff = cs_haz*Sminus
        return(riskRegression::rowCumSum(abs_risk_diff))
    })
    return(abs_risk)
}
## util for rewrapping list of lists 
rewrap <- function(ll, comb_fun = rbind){
    elements_in_list = length(ll[[1]])
    names_of_elements = names(ll[[1]])
    if(class(comb_fun) == "function")
        rep(c(comb_fun), times = elements_in_list)
    stopifnot(elements_in_list == length(comb_fun))
    out = lapply(1:elements_in_list, function(ii){
        do.call(comb_fun[[ii]], lapply(ll, function(xx){
            xx[[ii]]
        }))
    })
    names(out) = names_of_elements
    return(out)
}
## state learner
statelearner <- function(learners,
                         data,
                         time,
                         integrate = TRUE,
                         split.method = "cv5",
                         collapse = TRUE,
                         B=1,
                         verbose = FALSE,
                         time_grid_length = 100,
                         time_name = "time",
                         status_name = "status",
                         cause_codes = c("1" = 1, "2" = 2, "c" = 0),
                         vars = NULL){
    requireNamespace(package = "data.table")
    comp_event_present = !is.null(learners$cause2)
    ## NB: Expecting that time variable is called time
    partition = riskRegression::getSplitMethod(split.method=split.method,B=B,N=data[, .N])
    wd = form_data(data, time = time_name, status = status_name, cause_codes = cause_codes, vars = vars)
    setDT(wd)
    stopifnot(length(time) == 1)
    time_grid = seq(0, time, length.out = time_grid_length)
    ## fix names
    if(is.null(names(learners$cause1)))
        names(learners$cause1) = paste0("cause1_", 1:length(learners$cause1))
    if(is.null(names(learners$censor)))
        names(learners$censor) = paste0("censor_", 1:length(learners$censor))
    if(comp_event_present){
        if(is.null(names(learners$cause2)))
            names(learners$cause2) = paste0("cause2_", 1:length(learners$cause2))
    }
    list_cv_fits = lapply(1:partition$B, function(bb){
        if(verbose) message("--- Running on split ", bb, " out of ", B, " ---")
        raw_out0 = lapply(1:partition$k, function(fold_k){
            if(verbose) message("------ Running on split ", fold_k, " out of ", partition$k)
            train = wd[partition$index(bb)!=fold_k]
            test = wd[partition$index(bb)==fold_k]
            ## Construct time grid to evaluate hazard function
            eval_times = sort(unique(c(0, train[, sort(unique(time))])))
            eval_times = eval_times[eval_times <= max(time_grid)]
            ## Construct list of cumulative hazard predictions
            ## This might give memory problems...
            ## browser()
            time_fit = Sys.time()
            super_list_ch = lapply(names(learners), function(l_name){
                lapply(learners[[l_name]], function(mm){
                    fit = do.call(fit_cause_model,c(mm, list(data = train, cause = l_name)))
                    ch = predictCHF(fit, newdata = test, times = eval_times)
                    return(ch)
                })
            })
            ## ## Test
            ## fit_cause_model("GLMnet", train, "cause1")
            ## ## Test end            
            names(super_list_ch) = names(learners)
            list_ch_cause1 = super_list_ch[["cause1"]]
            list_ch_censor = super_list_ch[["censor"]]
            if(comp_event_present)
                list_ch_cause2 = super_list_ch[["cause2"]]
            time_fit = Sys.time() - time_fit
            ## Calculating Brier scores in holds out samples
            ## Construct matrices of events
            time_brier_calc = Sys.time()
            grid_mat = matrix(time_grid,ncol = length(time_grid),nrow = test[, .N],byrow = TRUE)
            event_mat = apply(grid_mat, 2, function(x) 1*(test[, time] <= x))
            cause1_event = apply(event_mat, 2, function(x) x*(test[, cause1 == 1]))
            censor_event = apply(event_mat, 2, function(x) x*(test[, censor == 1]))
            if(comp_event_present)
                cause2_event = apply(event_mat, 2, function(x) x*(test[, cause2 == 1]))
            if(comp_event_present)
                model_grid = expand.grid(cause1 = names((list_ch_cause1)),
                                         cause2 = names((list_ch_cause2)),
                                         censor = names((list_ch_censor)))
            else
                model_grid = expand.grid(cause1 = names((list_ch_cause1)),
                                         censor = names((list_ch_censor)))
            setDT(model_grid)
            brier_scores = do.call(rbind, lapply(1:nrow(model_grid), function(ii){
                base_dt = model_grid[ii]
                if(comp_event_present){
                    abs_risks = abs_risk_from_cschf(list_ch_cause1[[base_dt[, cause1]]],
                                                    list_ch_cause2[[base_dt[, cause2]]],
                                                    list_ch_censor[[base_dt[, censor]]])
                    event_counts = list(cause1_event, cause2_event, censor_event)
                }else{
                    abs_risks = abs_risk_from_cschf(list_ch_cause1[[base_dt[, cause1]]],
                                                    list_ch_censor[[base_dt[, censor]]])
                    event_counts = list(cause1_event, censor_event)
                }
                loss = sum(as.numeric(sapply(1:length(abs_risks), function(ll){
                    risk_pred = abs_risks[[ll]][, sindex(jump.times = eval_times,eval.times = time_grid)]
                    obs_event = event_counts[[ll]]
                    ## Check this ok!
                    brier = (risk_pred-obs_event)^2
                    if (integrate)
                        loss0 = mean(rowSums(brier*diff(time_grid)[1]))
                    else
                        loss0 = mean(brier[, dim(brier)[2]])
                    return(loss0)
                })))
                out = copy(base_dt)
                out[, loss := loss]
                return(out)                
            }))
            time_brier_calc = Sys.time() - time_brier_calc
            return(list(brier_scores = brier_scores,
                        time_fit =  time_fit,
                        time_brier_calc = time_brier_calc))
        })
        out0 = rewrap(raw_out0, c(rbind, sum, sum))
        if(comp_event_present)
            out0$brier_scores = out0$brier_scores[, .(loss=mean(loss)), .(cause1, cause2, censor)]
        else
            out0$brier_scores = out0$brier_scores[, .(loss=mean(loss)), .(cause1, censor)]
        out0$brier_scores[, b:=bb]
        return(out0)
    })
    rw_list_cv_fits = rewrap(list_cv_fits, c(rbind, sum, sum))
    cv_fits = rw_list_cv_fits$brier_scores
    if(verbose){
        all_time_fit = rw_list_cv_fits$time_fit
        all_time_brier_calc = rw_list_cv_fits$time_brier_calc
        message("Spend ", round(as.numeric(all_time_fit)/60, digits = 2), " minutes fitting models\n",
                "and   ", round(as.numeric(all_time_brier_calc)/60, digits = 2), " minutes calculating Brier scores.")
    }        
    if (B>1){
        if(comp_event_present)
            ave_cv_fit = cv_fits[, .(loss=mean(loss), sd = sd(loss)), .(cause1, cause2, censor)]
        else
            ave_cv_fit = cv_fits[, .(loss=mean(loss), sd = sd(loss)), .(cause1, censor)]
    }else{
        ave_cv_fit = cv_fits
    }
    if(collapse == TRUE)
        cv_return <- ave_cv_fit
    else
        cv_return <- cv_fits
    setkey(cv_return,loss)
    return(cv_return)
}


######################################################################
### statelearner-v2.R ends here
