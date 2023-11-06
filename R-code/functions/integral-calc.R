### integral-calc.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Oct 18 2023 (15:10) 
## Version: 
## Last-Updated: Nov  3 2023 (10:43) 
##           By: Anders Munch
##     Update #: 525
#----------------------------------------------------------------------
## 
### Commentary:
##
## These functions calculate varioues integrals
## based on piece-wise constant functions.
## 
### Change Log:
#----------------------------------------------------------------------
##

## Notes:
## Later on, consider implementing the integral functions where we supply matrices.
## Then they can be reused later on.
## Should minimize the number of times we need to call predict function.
## But then we need to export things in the end.
## Could we somehow write a completely functional expression, so that R figures out
## when to evaluate the supplied functions. 
## Thus only provide functions once, and number of chuncks once!
## All other functions should the be symboblic combinations of the supplied functions.
## Would also make stuff easier to read.
## And then everything the combiner (chunk) function could be exported to C++.

### Code:
library(data.table)
library(prodlim)
library(riskRegression)

## Meta functions:
minus_t_fun <- function(Lambda, jump_points){
    t_minus = min(diff(sort(unique(jump_points))))/2
    fun0 = function(newdata, times){
        times0 = pmax(0, times-t_minus)
        Lambda(newdata, times0)
    }
    return(fun0)
}
construct_pred_fun <- function(object, ...){
    UseMethod("construct_pred_fun",object)
}
construct_pred_fun.coxph <- function(model, ...){
    out = function(newdata,times)
        matrix(predictCox(model, newdata = newdata, times = times, type = "cumhazard")$cumhazard, ncol = length(times))
    return(out)
}
construct_pred_fun.glm <- function(model, ...){
    out = function(newdata)
        matrix(predict(model, newdata = newdata, type = "response"), ncol = 1)
    return(out)
    
}
#### at risk funcion
at_risk_fun <- function(newdata,times,time_name = "time"){
    spring_times = newdata[[time_name]]
    grid_mat = matrix(times,ncol = length(times),nrow = length(spring_times),byrow = TRUE)
    ## event_mat = apply(grid_mat, 2, function(x) 1*(spring_times >= x))
    event_mat = do.call(cbind, lapply(1:ncol(grid_mat), function(xx) 1*(spring_times >= grid_mat[, xx])))
    return(event_mat)
}


#### Helper functions
## Format input and send messages to user
intern_setup_message <- function(start,end,jump_points,chunks,integrate,verbose = TRUE){
    stopifnot(start <= end)
    jump_points = sort(unique(c(start, jump_points, end)))
    if(!(length(jump_points)>1))
        return(list(status = "stop"))
    wtimes = jump_points[start <= jump_points & jump_points <= end]
    n_t = length(wtimes)
    if(chunks >= n_t){
        chunks = (n_t-1)
        if(verbose)
            message(paste0("Number of chuncks is bigger than number of relevant jump_points. Chunks set to ", n_t-1, "."))
    }
    if(chunks>1){
        if(!integrate){
            if(verbose)
                message("Need to integrate when chunks > 1, so the 'integrate' argument is ignored.")
            integrate = TRUE
        }
    }
    return(list(status = "continue", wtimes = wtimes, chunks = chunks, integrate = integrate))
}
## Calculate integral, assuming correctly formatted time points.
intern_int <- function(A,
                       B,
                       start,
                       end,
                       data,
                       wtimes,
                       integrate,
                       report_raw){
    ## Calculate the integral/sum:
    ## \int_start^end A(s) B(ds)
    ##
    A0 = A(newdata = data,times = wtimes[-1])
    B0 = B(newdata = data,times = wtimes)
    B0_diff = do.call(rbind, lapply(1:nrow(B0), function(xx) diff(B0[xx,])))
    if(report_raw)
        return(list(A0 = A0, B0_diff = B0_diff))
    integrands = A0 * B0_diff
    if(integrate)
        out = rowSums(integrands)
    else
        out = integrands
    return(out)
}
## Running integral, assuming correctly formatted time points.
intern_running_int <- function(A,B,start,end,jump_points){
    wtimes = intern_setup_message(start = start, end = end, jump_points = jump_points, chunks = 1, integrate = FALSE, verbose = FALSE)$wtime
    out_fun = function(newdata, times){
        integrands = intern_int(A = A,
                                B = B,
                                start = start,
                                end = end,
                                wtimes = wtimes,
                                data = newdata,
                                integrate = FALSE,
                                report_raw = FALSE)
        fun_grid = do.call(rbind, lapply(1:nrow(integrands), function(xx) cumsum(integrands[xx,])))
        ## Must be 0 at time 0, because it is \int_start^start ...
        fun_grid = cbind(0, fun_grid)        
        out = fun_grid[, prodlim::sindex(jump.times = wtimes, eval.times = times)]
        return(out)
    }
    return(out_fun)
}

#### Integrals
## Divides calculation into number of chunks to save memory.
integral1 <- function(A,
                      B,
                      start = 0,
                      end,
                      data,
                      jump_points,
                      chunks = 1,
                      integrate = TRUE){
    ## Calculate the integral/sum:
    ## \int_start^end A(s) B(ds)
    setup0 = intern_setup_message(start = start, end = end,jump_points = jump_points,chunks = chunks, integrate = integrate)
    if(setup0$status == "stop")
        return(numeric(nrow(data))) ## Just return a vector of zeros
    wtimes = setup0$wtimes
    chunks = setup0$chunks
    integrate = setup0$integrate
    if(chunks == 1){
        out = intern_int(A = A,B = B,start = start,end = end,data = data,wtimes = wtimes,integrate = integrate, report_raw = FALSE)
    }else{
        ## chop up integral to save memory:
        chop_points = quantile(wtimes, seq(0,1,length.out = chunks+1), type = 1)
        out = numeric(nrow(data))
        for(ii in 1:(length(chop_points)-1)){
            start_ii = chop_points[ii]
            end_ii = chop_points[ii+1]
            times_ii = wtimes[start_ii <= wtimes & wtimes <= end_ii]
            int_ii = intern_int(A = A,B = B,start = start_ii,end = end_ii,data = data,wtimes = times_ii,integrate = TRUE, report_raw = FALSE)
            out = out + int_ii
        }
    }
    return(out)
}
integral2 <- function(A, B, C, D, start = 0, end, data, jump_points, chunks = 1, integrate = TRUE){
    ## Calculate the double integral/sum:
    ## \int_start^end C(t) [\int_0^t A(s) B(ds)] D(dt)
    setup0 = intern_setup_message(start = start, end = end,jump_points = jump_points,chunks = chunks, integrate = integrate)
    if(setup0$status == "stop")
        return(numeric(nrow(data))) ## Just return a vector of zeros
    wtimes = setup0$wtimes
    chunks = setup0$chunks
    integrate = setup0$integrate
    if(chunks == 1){
        inner_int_fun = intern_running_int(A = A, B = B, start = start, end = end, jump_points = wtimes)
        A_outer = function(newdata,times) (inner_int_fun(newdata,times) * C(newdata,times))
        out = intern_int(A = A_outer,B = D,start = start,end = end,data = data,wtimes = wtimes,integrate = integrate,report_raw = FALSE)
    }else{
        chop_points = quantile(wtimes, seq(0,1,length.out = chunks+1), type = 1)
        out = numeric(nrow(data))
        inner_int_last = numeric(nrow(data))
        for(ii in 1:(length(chop_points)-1)){
            start_ii = chop_points[ii]
            end_ii = chop_points[ii+1]
            times_ii = wtimes[start_ii <= wtimes & wtimes <= end_ii]
            inner_int_fun = intern_running_int(A = A, B = B, start = start_ii, end = end_ii, jump_points = times_ii)
            ## A_outer = function(newdata, times) {C(newdata, times)*(inner_int_fun(newdata,times) + inner_int_last)}
            A_outer = function(newdata, times) {inner_int_fun(newdata,times) + inner_int_last}
            outer_int_raw = intern_int(A = A_outer,B = D,start = start_ii,end = end_ii,data = data,wtimes = times_ii,integrate = FALSE,report_raw = TRUE)
            ## Manual hack
            C0 = C(newdata = data, times = times_ii[-1])
            out = out + rowSums(C0 * outer_int_raw$A0 * outer_int_raw$B0_diff)
            ## Update inner integral up to this time point
            inner_int_last = outer_int_raw$A0[,ncol(outer_int_raw$A0)]
        }
    }
    return(out)
}
integral0_itime <- function(A,
                            times,
                            data,
                            chunks = 1){
    n_t = length(times)
    stopifnot(n_t == nrow(data))
    chop_ii = seq(1,n_t,by = max(5, ceiling(n_t/(chunks+1))))
    if((length(chop_ii) <= 2) | (chunks == 1))
        return(diag(A(newdata = data, times = times)))
    chop_ii[length(chop_ii)] = n_t+1
    db_times = data.table(id_orig = 1:n_t, time = times)
    setorder(db_times, time)
    db_times[, id_sorted := 1:n_t]    
    out_sorted = do.call(c, lapply(1:(length(chop_ii)-1), function(ii){
        ii_start = chop_ii[ii]
        ii_end = chop_ii[ii+1]-1
        times_ii = db_times[ii_start:ii_end, time]
        data_ii = data[db_times[ii_start:ii_end, id_orig]]
        return(diag(A(newdata = data_ii, times = times_ii)))
    }))
    setorder(db_times, id_orig)
    out = out_sorted[db_times[, id_sorted]]
    return(out)
}
## Construct integral for
## \int_start^end A(s) N(\diff s),
## where N is a counting process.
count_integral <- function(A,
                           t,
                           times,
                           events,
                           data,
                           chunks = 1){
    stopifnot((length(times) == nrow(data)) & (length(events) == nrow(data)))
    out = numeric(nrow(data))
    rel_ind = as.logical((events == 1) & (times <= t))
    if(sum(rel_ind) != 0){
        wd = copy(data[rel_ind])
        wtimes = times[rel_ind]
        out[rel_ind] = integral0_itime(A = A,times = wtimes,data = wd,chunks = chunks)
    }
    return(out)
}


#### Old versions
## internal_integral1 <- function(A,
##                                 B,
##                                 end,
##                                 start = 0,
##                                 data,
##                                 jump_points,
##                                 integrate = TRUE){
##     ## Calculate the integral/sum:
##     ## \int_start^end A(s) B(ds)
##     ##
##     jump_points = sort(unique(c(start, jump_points, end)))
##     stopifnot(length(jump_points)>1)
##     wtimes = jump_points[start <= jump_points & jump_points <= end]
##     n_t = length(wtimes)
##     A0 = A(newdata = data,times = wtimes)
##     B0 = B(newdata = data,times = wtimes)
##     ## B0_diff = t(apply(B0, 1, diff))
##     B0_diff = do.call(rbind, lapply(1:nrow(B0), function(xx) diff(B0[xx,])))    
##     ## Why not:
##     integrands = A0[, -1] * B0_diff
##     if(integrate)
##         out = rowSums(integrands)
##     else
##         out = integrands
##     return(out)
## }

## internal_integral2 <- function(A, B, C, D, end, start = 0, data, jump_points, integrate = TRUE){
##     ## Calculate the double integral/sum:
##     ## \int_start^end C(t) [\int_0^t A(s) B(ds)] D(dt)
##     ##
##     jump_points = sort(unique(c(start, jump_points, end)))
##     wtimes = jump_points[start <= jump_points & jump_points <= end]    
##     n_t = length(wtimes)
##     mf0 = function(newdata, times){
##         integrands = internal_intergral1_v2(A = A, B = B, start = start, end = end, jump_points = times, newdata = newdata, integrate = FALSE)
##         out = do.call(rbind, lapply(1:nrow(integrands), function(xx) cumsum(integrands[xx,])))
##         ## Must be 0 at time 0, because it is \int_start^start ...
##         out = cbind(0, out)
##         return(out)
##     }
##     A_meta = function(newdata,times) (mf0(newdata,times) * C(newdata,times))
##     outer_int = internal_intergral1_v2(A = A_meta,
##                                        B = D,
##                                        end = end,
##                                        start = start,
##                                        newdata = data,
##                                        jump_points = wtimes,
##                                        integrate = integrate)
##     return(outer_int)
## }


## internal_intergral1 <- function(A,
##                                 B = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)),
##                                 C,
##                                 end,
##                                 start = 0,
##                                 newdata,
##                                 jump_points,
##                                 integrate = TRUE){
##     ## Calculate the integral/sum:
##     ## \int_start^end exp(A(s)) * B(s) C(ds)
##     ##
##     jump_points = sort(unique(c(start, jump_points, end)))
##     stopifnot(length(jump_points)>1)
##     wtimes = jump_points[start <= jump_points & jump_points <= end]
##     n_t = length(wtimes)
##     A0 = A(newdata = newdata,times = wtimes)
##     B0 = B(newdata = newdata,times = wtimes)
##     C0 = C(newdata = newdata,times = wtimes)
##     C0_diff = t(apply(C0, 1, diff))
##     ## integrands = exp(A0[, -length(n_t)]) * B0[, -length(n_t)] * C0_diff
##     ## Why not:
##     integrands = exp(A0[, -1]) * B0[, -1] * C0_diff
##     if(integrate)
##         out = rowSums(integrands)
##     else
##         out = integrands
##     return(out)
## }

## internal_intergral2 <- function(A = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)), B, C, D, end, start = 0, newdata, jump_points, integrate = TRUE){
##     ## Calculate the double integral/sum:
##     ## \int_start^end A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
##     ##
##     jump_points = sort(unique(c(start, jump_points, end)))
##     wtimes = jump_points[start <= jump_points & jump_points <= end]    
##     n_t = length(wtimes)
##     mf0 = function(newdata, times){
##         integrands = internal_intergral1(A = B, C = C, start = start, end = end, jump_points = times, newdata = newdata, integrate = FALSE)
##         out = do.call(rbind, lapply(1:nrow(integrands), function(xx) cumsum(integrands[xx,])))
##         ## Must be 0 at time 0, because it is \int_start^start ...
##         out = cbind(0, out)
##         return(out)
##     }
##     outer_int = internal_intergral1(A = function(newdata,times) matrix(0,nrow = nrow(newdata),ncol = length(times)),
##                                     B = mf0,
##                                     C = D,
##                                     end = end,
##                                     start = start,
##                                     newdata = newdata,
##                                     jump_points = wtimes,
##                                     integrate = integrate)
##     return(outer_int)
## }


## internal_intergral2_check <- function(A = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)),
##                                       B,
##                                       C,
##                                       D,
##                                       end,
##                                       start = 0,
##                                       newdata,
##                                       jump_points,
##                                       integrate = TRUE){
##     ## Calculate the double integral/sum:
##     ## \int_start^end A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
##     ##
##     jump_points = sort(unique(c(start, jump_points, end)))
##     wtimes = jump_points[start <= jump_points & jump_points <= end]    
##     n_t = length(wtimes)
##     B0 = B(newdata = newdata,times = wtimes)
##     C0 = C(newdata = newdata,times = wtimes)
##     C0_diff = t(apply(C0, 1, diff))
##     meta_B0_int = exp(B0[, -1]) * C0_diff
##     meta_B0 = do.call(rbind, lapply(1:nrow(meta_B0_int), function(xx) cumsum(c(meta_B0_int[xx,]))))
##     A0 = A(newdata = newdata,times = wtimes)
##     D0 = D(newdata = newdata,times = wtimes)
##     D0_diff = t(apply(D0, 1, diff))
##     integrands = A0[, -1] * meta_B0 * D0_diff
##     if(integrate)
##         out = rowSums(integrands)
##     else
##         out = integrands
##     return(out)
## }


######################################################################
### integral-calc.R ends here
