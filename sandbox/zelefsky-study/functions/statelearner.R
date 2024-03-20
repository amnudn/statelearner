### statelearner.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2023 (13:29) 
## Version: 
## Last-Updated: Aug 21 2023 (10:12) 
##           By: Anders Munch
##     Update #: 89
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
statelearner <- function(learners,
                         data,
                         times,
                         integrate = TRUE,
                         split.method = "cv5",
                         collapse = TRUE,
                         B=1,
                         verbose = FALSE,
                         time_grid_length = 100){
    requireNamespace(package = "data.table")
    ## NB: Expecting that time variable is called time
    partition = riskRegression::getSplitMethod(split.method=split.method,B=B,N=data[, .N])
    wd = data.table::copy(data)
    setDT(wd)
    wd[, event_c := 1*(!event)]
    stopifnot(length(times)>0)
    if (length(times) == 1)
        time_grid = seq(0, times, length.out = time_grid_length)
    else
        time_grid = times
    ## fix names
    if(is.null(names(learners$state)))
        names(learners$state) = paste0("out_", 1:length(learners$state))
    if(is.null(names(learners$censoring)))
        names(learners$censoring) = paste0("cens_", 1:length(learners$censoring))
    cv_fits = do.call(rbind, lapply(1:partition$B, function(bb){
        if(verbose) message("--- Running on split ", bb, " out of ", B, " ---")
        out0 = do.call(rbind, lapply(1:partition$k, function(fold_k){
            if(verbose) message("------ Running on split ", fold_k, " out of ", partition$k)
            train = wd[partition$index(bb)!=fold_k]
            test = wd[partition$index(bb)==fold_k]
            ## Construct time grid to evaluate hazard function
            eval_times = sort(unique(c(0, train[, sort(unique(time))])))
            eval_times = eval_times[eval_times <= max(time_grid)]
            ## Construct list of cumulative hazard predictions
            ## This might give memory problems...
            ## Solving this partly by only calculating all CHs for censoring only
            ## The other calculation happens in outer loop below
            ## -- then we have only L+1 matrices of CHs, instead of 2*L matrices
            list_ch_cens = lapply(learners$censoring, function(mm){
                fit = update(mm, data = train)
                ch = -log(1-predictRisk(fit, newdata = test, times = eval_times, product.limit = 0))
                return(ch)
            })
            ## Calculating Brier scores in holds out samples
            ## Construct matrices of events
            grid_mat = matrix(time_grid,ncol = length(time_grid),nrow = test[, .N],byrow = TRUE)
            event_mat = apply(grid_mat, 2, function(x) 1*(test[, time] <= x))
            status_out = apply(event_mat, 2, function(x) x*(test[, event]))
            status_cens = apply(event_mat, 2, function(x) x*(test[, event_c]))
            train_times = train[, sort(unique(time))]
            ## browser()
            brier_scores = do.call(rbind, lapply(names(learners$state), function(out_n){
                fit_out = update(learners$state[[out_n]], data = train)
                ch_out = -log(1-predictRisk(fit_out, newdata = test, times = eval_times, product.limit = 0))
                do.call(rbind, lapply(names(learners$censoring), function(cens_n){
                    ch_cens = list_ch_cens[[cens_n]]
                    ## abs_risk_from_ch <- function(ch1, ch2, time_grid, fit_times){
                    # survival just before time t
                    S = exp(-ch_out-ch_cens)
                    Sminus = cbind(1,S[,-ncol(S)])
                    h_out = t(apply(cbind(0,ch_out), 1, diff))
                    h_cens = t(apply(cbind(0,ch_cens), 1, diff))
                    ## NB, these are not abs risk, but the derivatives of them!
                    abs_risk_out = Sminus*h_out
                    abs_risk_cens = Sminus*h_cens
                    # integrating over time 
                    report_index = prodlim::sindex(jump.times = train_times,eval.times = time_grid)
                    ## Bad practice to use same name here.
                    abs_risk_out = riskRegression::rowCumSum(abs_risk_out[,1:max(report_index)])
                    abs_risk_cens = riskRegression::rowCumSum(abs_risk_cens[,1:max(report_index)])
                    ## if (any(is.na(abs_risk_out)) | any(is.na(abs_risk_cens))) browser() #stop("Missing values in predicted risks for censoring learner: ",cens_n)
                    fix_NaNs = function(ar_check, ar_replace){

                        fixme = which(apply(ar_check,1,function(u)sum(is.na(u)|is.infinite(u)))>0)
                        if (length(fixme)>0){ 
                            for (i in fixme){
                                tmp = ar_check[i,]
                                fix_col = (is.na(tmp)|is.infinite(tmp))
                                ar_check[i,fix_col] = ar_replace[i,fix_col]
                            }
                        }
                        return(ar_check)
                    }
                    abs_risk_cens = fix_NaNs(abs_risk_cens, 1-abs_risk_out)
                    abs_risk_out = fix_NaNs(abs_risk_out, 1-abs_risk_cens)
                    res1 = (cbind(0,abs_risk_out)[,1+report_index]-status_out)
                    res2 = (cbind(0,abs_risk_cens)[,1+report_index]-status_cens)
                    brier = 2*(res1^2 + res2^2 +(res1*res2))
                    ## if (any(is.na(brier))) browser() #stop("Missing values in predicted risks for censoring learner: ",cens_n)
                    if (integrate)
                        return(data.table(out_model = out_n,
                                          cens_model = cens_n,
                                          brier_score = mean(brier)))
                    else
                        return(data.table(out_model = out_n,
                                          cens_model = cens_n,
                                          brier_score = mean(brier[, dim(brier)[2]])))
                }))
            }))
            return(brier_scores)
        }))
        out0 = out0[, .(loss=mean(brier_score)), .(out_model, cens_model)]
        out0[, b:=bb]
        return(out0)
    }))
    if (B>1)
        ave_cv_fit = cv_fits[, .(loss=mean(loss), sd = sd(loss)), .(out_model, cens_model)]
    else
        ave_cv_fit = cv_fits
    if(collapse == TRUE)
        cv_return <- ave_cv_fit
    else
        cv_return <- cv_fits
    setkey(cv_return,loss)
    return(cv_return)
}


######################################################################
### statelearner.R ends here
