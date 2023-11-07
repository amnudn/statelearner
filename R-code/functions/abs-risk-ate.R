### martingale-int.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Oct 27 2023 (16:08) 
## Version: 
## Last-Updated: Nov  7 2023 (11:22) 
##           By: Anders Munch
##     Update #: 426
#----------------------------------------------------------------------
## 
### Commentary:
##
## Implements a one-step estimator of the average treatment effect
## on cause 1.
##
## Calculation is based on fitted function that calculate cumulative hazard function.
##
## See implementation doc for details on the nomenclature.
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(here)
library(data.table)

## naiv (data, t, Lambda1, Lambda2, jump_points)
naiv <- function(data, t, Lambda1, Lambda2, jump_points, chunks = 1){
    Lambda1_ = minus_t_fun(Lambda1,jump_points)
    Lambda2_ = minus_t_fun(Lambda2,jump_points)
    A_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
    data_treat0 = copy(data)[, A := 0]
    data_treat1 = copy(data)[, A := 1]
    preds0 = integral1(A = A_fun,B = Lambda1,data = data_treat0,end = t,jump_points = jump_points, chunks = chunks)
    preds1 = integral1(A = A_fun,B = Lambda1,data = data_treat1,end = t,jump_points = jump_points, chunks = chunks)
    out = list("A=1" = preds1, "A=0" = preds0)
    return(out)
}
termW <- function(data, pi){
    pi_fit = as.numeric(pi(data))
    out1 = numeric(nrow(data))
    A1_ind = data[, A == 1]
    out1[A1_ind] = 1/pi_fit[A1_ind]    
    out0 = numeric(nrow(data))
    out0[!A1_ind] = 1/(1-pi_fit[!A1_ind])
    out = list("A=1" = out1, "A=0" = out0)
    return(out)
}
termA <- function(data, t, Lambda1, Gamma, jump_points,chunks = 1){
    Gamma_ = minus_t_fun(Gamma,jump_points)
    ## Counting term
    out_count = count_integral(A = function(newdata,times) exp(Gamma_(newdata,times)),t = t,times = data[, time],events = data[, (status == 1)],data = data,chunks = chunks)
    ## Compensator
    A_fun = function(newdata,times) {exp(Gamma_(newdata,times))*at_risk_fun(newdata, times)}
    out_comp = integral1(A = A_fun,B = Lambda1,end = t,data = data,jump_points = jump_points,chunks = chunks)
    out = out_count-out_comp
    return(out)
}
termB <- function(data, t, Lambda1, Lambda2, Gamma, jump_points,chunks = 1){
    Lambda1_ = minus_t_fun(Lambda1,jump_points)
    Lambda2_ = minus_t_fun(Lambda2,jump_points)
    Gamma_ = minus_t_fun(Gamma,jump_points)
    g_term = integral1(A = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times))),
                       B = Lambda1,
                       end = t,
                       data = data,
                       jump_points = jump_points,
                       chunks = chunks)
    ## Counting term
    out_count = count_integral(A = function(newdata,times) exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times)),
                               t = t,                                
                               times = data[, time],
                               events = data[, (status>0)],
                               data = data,
                               chunks = chunks)
    ## Compensator term
    A_fun = function(newdata,times) {
        exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
            at_risk_fun(newdata, times)
    }
    out_comp = integral1(A = A_fun,
                         B = function(newdata,times) {Lambda1(newdata,times)+Lambda2(newdata,times)},
                         end = t,
                         data = data,
                         jump_points = jump_points,
                         chunks = chunks)
    out = g_term*(out_count-out_comp)
    return(out)
}
termC <- function(data, t, Lambda1, Lambda2, Gamma, jump_points,chunks = 1){
    Lambda1_ = minus_t_fun(Lambda1,jump_points)
    Lambda2_ = minus_t_fun(Lambda2,jump_points)
    Gamma_ = minus_t_fun(Gamma,jump_points)
    ## Counting term
    out_count = count_integral(A = function(newdata,times)  exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times)),
                               t = t,                                
                               times = data[, time],
                               events = data[, (status>0)],
                               data = data,
                               chunks = chunks)
    ## Only calculate following for non-zero terms to minimize positivite issues
    zero_count = (out_count == 0)
    if(!all(zero_count)){
        wd = copy(data[!zero_count])
        g_term = integral1(A = function(newdata,times) {at_risk_fun(newdata,times)*exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))},
                           B = Lambda1,
                           end = t,
                           data = wd,
                           jump_points = jump_points,
                           chunks = chunks)
        out_count[!zero_count] = out_count[!zero_count]*g_term
    }
    ## Compensator term
    C_fun = function(newdata,times) {
        exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
            at_risk_fun(newdata, times)
    }
    out_comp = integral2(A = function(newdata,times){exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))},
                         B = Lambda1,
                         C = C_fun,
                         D = function(newdata,times) {Lambda1(newdata,times)+Lambda2(newdata,times)},
                         end = t,
                         data = data,
                         jump_points = jump_points,
                         chunks = chunks)
    out = (out_count-out_comp)
    return(out)
}
raw_os_abs_risk_ate <- function(data, t, Lambda1, Lambda2, Gamma, pi, jump_points,chunks = 1,collapse = TRUE){
    W_i = termW(data = data, pi = pi)
    A_i = termA(data = data, t = t, Lambda1 = Lambda1, Gamma = Gamma, jump_points = jump_points, chunks = chunks)
    B_i = termB(data = data, t = t, Lambda1 = Lambda1, Lambda2 = Lambda2, Gamma = Gamma, jump_points = jump_points, chunks = chunks)
    C_i = termC(data = data, t = t, Lambda1 = Lambda1, Lambda2 = Lambda2, Gamma = Gamma, jump_points = jump_points, chunks = chunks)
    naiv_i = naiv(data = data, t = t, Lambda1 = Lambda1, Lambda2 = Lambda2, jump_points = jump_points, chunks = chunks)
    out = list(naiv1_i = naiv_i[["A=1"]],
               naiv0_i = naiv_i[["A=0"]],
               W1_i = W_i[["A=1"]],
               W0_i = W_i[["A=0"]],
               A_i = A_i, B_i = B_i, C_i = C_i)
    if(class(collapse) == "logical")
        collapse = 2*collapse
    if(collapse)
        out = list("A=1" = out$W1_i*(A_i - B_i + C_i) + out$naiv1_i,
                   "A=0" = out$W0_i*(A_i - B_i + C_i) + out$naiv0_i)
    if(collapse == 2){
        out = list("naiv" =
                       list("A=1" = mean(naiv_i[["A=1"]]),
                            "A=0" = mean(naiv_i[["A=0"]]),
                            "ATE" = mean(naiv_i[["A=1"]]-naiv_i[["A=0"]])),
                   "one_step" =
                       list("A=1" = mean(out[["A=1"]]),
                            "A=0" = mean(out[["A=0"]]),
                            "ATE" = mean(out[["A=1"]])-mean(out[["A=0"]])))
    }
    return(out)
}
os_abs_risk_ate <- function(data, eval_times, fit_1, fit_2, fit_cens, fit_treat, jump_points = data[, time],chunks = 1){
    L1 = function(newdata,times) predictCHF(fit_1, newdata, times)
    L2 = function(newdata,times) predictCHF(fit_2, newdata, times)
    G = function(newdata,times) predictCHF(fit_cens, newdata, times)
    pi = function(newdata) predictTreat(fit_treat, newdata)
    cause_est = list(cause1 = L1, cause2 = L2)
    out = do.call(rbind, lapply(eval_times, function(tt){
        do.call(rbind, lapply(1:length(cause_est), function(ii){
            cause_interest = names(cause_est)[ii]
            L1_ii = cause_est[[ii]]
            L2_ii = cause_est[[1+(ii %% 2)]]
            raw0 = raw_os_abs_risk_ate(data = data,
                                       t = tt,
                                       Lambda1 = L1_ii,
                                       Lambda2 = L2_ii,
                                       Gamma = G,
                                       pi = pi,
                                       jump_points = jump_points,
                                       chunks = chunks,
                                       collapse = 0)
            naiv1 = mean(raw0$naiv1_i)
            naiv0 = mean(raw0$naiv0_i)
            debias_term1 = with(raw0, mean(W1_i*(A_i - B_i + C_i)))
            debias_term0 = with(raw0, mean(W0_i*(A_i - B_i + C_i)))
            see_1 = with(raw0, sd(naiv1_i + W1_i*(A_i - B_i + C_i))/sqrt(nrow(data)))
            see_0 = with(raw0, sd(naiv0_i + W0_i*(A_i - B_i + C_i))/sqrt(nrow(data)))
            see_ate = with(raw0, sd(naiv1_i + W1_i*(A_i - B_i + C_i) - (naiv0_i + W0_i*(A_i - B_i + C_i)))/sqrt(nrow(data)))
            effect_names = c("A=1","A=0","ATE")
            out0 = data.table(cause = cause_interest,
                              time = tt,
                              effect = effect_names,
                              est = c(naiv1+debias_term1, naiv0+debias_term0, naiv1+debias_term1 - (naiv0+debias_term0)),
                              see = c(see_1, see_0, see_ate))
        }))
    }))
    out[, ":="(lower = est-1.96*see, upper = est+1.96*see)]
    return(out[])
}


#### Old versions
## naiv_est <- function(data, t, Lambda1, Lambda2, jump_points, collapse = TRUE){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     A_fun = function(newdata,times) {-(Lambda1_(newdata,times) + Lambda2_(newdata,times))}
##     data_treat0 = copy(data)[, A := 0]
##     data_treat1 = copy(data)[, A := 1]
##     preds0 = internal_intergral1(A = A_fun,C = Lambda1,newdata = data_treat0,end = t,jump_points = jump_points)
##     preds1 = internal_intergral1(A = A_fun,C = Lambda1,newdata = data_treat1,end = t,jump_points = jump_points)
##     if(collapse)
##         return(mean(preds1 - preds0))
##     else
##         return(data.table(Y1_init_est = preds1, Y1_init_est = preds0))
## }

## Q1_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE, chunks = 1){
##     Q_out = numeric(nrow(data))
##     rel_ind = data[,(status>0) & (time <= t)]
##     if(sum(rel_ind) == 0)
##         return(Q_out)
##     wd = copy(data[rel_ind])
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)
##     A_fun = function(newdata,times) {-(Lambda1_(newdata,times) + Lambda2_(newdata,times))}
##     not_at_risk0 = function(newdata, times) {1-at_risk_fun(newdata, times)}
##     term_w = exp(diag(Lambda1_(wd,wd[,time]) + Lambda2_(wd,wd[,time]) + Gamma_(wd,wd[,time])))
##     term_int = internal_intergral1(A = A_fun, B = not_at_risk0, C = Lambda1, end = t, newdata = wd, jump_points = jump_points)
##     Q_out[rel_ind] =  term_w*term_int
##     if(collapse)
##         return(mean(Q_out))
##     else
##         return(Q_out)
## }

##
## Q2_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)    
##     C_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
##     A_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     B_fun = function(newdata,times) {Lambda1(newdata,times) + Lambda2(newdata,times)}
##     out = integral2(A = A_fun,
##                     B = B_fun,
##                     C = C_fun,
##                     D = Lambda1,
##                     end = t,
##                     data = data,
##                     jump_points = jump_points,
##                     chunks = chunks)
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }


## ## Q1
## Q1_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE, chunks = 1, check = FALSE){
##     out = numeric(nrow(data))
##     rel_ind = data[,(status>0) & (time <= t)]
##     if(sum(rel_ind) == 0)
##         return(out)
##     wd = copy(data[rel_ind])
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)
##     term_w = integral0_itime(A = function(newdata,times) exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times)),
##                              times = wd[,time],
##                              data = wd,
##                              chunks = chunks)
##     A_fun = function(newdata,times) {
##         exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))*
##             (1-at_risk_fun(newdata, times))
##     }
##     term_int = integral1(A = A_fun, B = Lambda1, end = t, data = wd, jump_points = jump_points, chunks = chunks)
##     out[rel_ind] =  term_w*term_int
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }

## ## Q2
## Q2_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)    
##     h_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
##     gamma_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     bar_fun = function(newdata,times) {Lambda1(newdata,times) + Lambda2(newdata,times)}
##     term_plus = integral1(A = h_fun,
##                           B = Lambda1,
##                           end = t,
##                           data = data,
##                           jump_points = jump_points,
##                           chunks = chunks)*
##         integral1(A = gamma_fun,
##                   B = bar_fun,
##                   end = t,
##                   data = data,
##                   jump_points = jump_points,
##                   chunks = chunks)
##     term_minus = integral2(A = h_fun,
##                            B = Lambda1,
##                            C = gamma_fun,
##                            D = bar_fun,
##                            end = t,
##                            data = data,
##                            jump_points = jump_points,
##                            chunks = chunks)
##     out = term_plus-term_minus
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }


## Q2_check <- function(data, t, Lambda1, Lambda2, Gamma, jump_points,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)    
##     C_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
##     A_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     ## B_fun = function(newdata,times) {Lambda1(newdata,times) + Lambda2(newdata,times)}
##     ## Test:
##     B_fun = function(newdata,times) {Lambda1_(newdata,times) + Lambda2_(newdata,times)}
##     out = integral2(A = A_fun,
##                     B = B_fun,
##                     C = C_fun,
##                     D = Lambda1,
##                     end = t,
##                     data = data,
##                     jump_points = jump_points,
##                     chunks = chunks,
##                     integrate = FALSE)
##     ## Alternative form:
##     h_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
##     gamma_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     bar_fun = function(newdata,times) {Lambda1(newdata,times) + Lambda2(newdata,times)}
##     term_plus1 = integral1(A = h_fun,B = Lambda1,end = t,data = data,jump_points = jump_points,chunks = chunks,integrate = FALSE)    
##     term_plus2 = integral1(A = gamma_fun,B = bar_fun,end = t,data = data,jump_points = jump_points,chunks = chunks,integrate = FALSE)
##     term_minus = integral2(A = h_fun,B = Lambda1,C = gamma_fun,D = bar_fun,end = t,data = data,jump_points = jump_points,chunks = chunks,integrate = FALSE)
##     return(list(out = out, term_plus1 = term_plus1, term_plus2 = term_plus2, term_minus = term_minus))
## }

## ## Checking:
## R1_count_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)
##     g_term = integral1(A = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times))),
##                        B = Lambda1,
##                        end = t,
##                        data = data,
##                        jump_points = jump_points,
##                        chunks = chunks)
##     ## Counting term
##     out_count = numeric(nrow(data))
##     rel_ind = data[,(status>0) & (time <= t)]
##     if(sum(rel_ind) != 0){
##         wd = copy(data[rel_ind])
##         out_count[rel_ind] = integral0_itime(A = function(newdata,times) exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times)),
##                                              times = wd[,time],
##                                              data = wd,
##                                              chunks = chunks)
##     }
##     out = g_term*out_count
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }
## R1_comp_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)
##     g_term = integral1(A = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times))),
##                        B = Lambda1,
##                        end = t,
##                        data = data,
##                        jump_points = jump_points,
##                        chunks = chunks)
##     ## Compensator term
##     A_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     out_comp = integral1(A = A_fun,
##                          B = function(newdata,times) {Lambda1(newdata,times)+Lambda2(newdata,times)},
##                          end = t,
##                          data = data,
##                          jump_points = jump_points,
##                          chunks = chunks)
##     out = g_term*out_comp
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }

## R2_count_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)
##     ## Counting term
##     out_count = numeric(nrow(data))
##     rel_ind = data[,(status>0) & (time <= t)]
##     if(sum(rel_ind) != 0){
##         wd = copy(data[rel_ind])
##         g_term = integral1(A = function(newdata,times) {at_risk_fun(newdata,times)*exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))},
##                            B = Lambda1,
##                            end = t,
##                            data = wd,
##                            jump_points = jump_points,
##                            chunks = chunks)
##         int_term = integral0_itime(A = function(newdata,times) exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times)),
##                                    times = wd[,time],
##                                    data = wd,
##                                    chunks = chunks)
##         out_count[rel_ind] = g_term*int_term
##     }
##     out = out_count
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }
## R2_comp_fun <- function(data, t, Lambda1, Lambda2, Gamma, jump_points, collapse = TRUE,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)
##     ## Compensator term
##     C_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     out_comp = integral2(A = function(newdata,times){exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))},
##                          B = Lambda1,
##                          C = C_fun,
##                          D = function(newdata,times) {Lambda1(newdata,times)+Lambda2(newdata,times)},
##                          end = t,
##                          data = data,
##                          jump_points = jump_points,
##                          chunks = chunks)
##     out = out_comp
##     if(collapse)
##         return(mean(out))
##     else
##         return(out)
## }

## ## pi_w (data, pi)



## function(data, t, Lambda1, Lambda2, Gamma, jump_points,chunks = 1){
##     Lambda1_ = minus_t_fun(Lambda1,jump_points)
##     Lambda2_ = minus_t_fun(Lambda2,jump_points)
##     Gamma_ = minus_t_fun(Gamma,jump_points)    
##     C_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
##     A_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     ## B_fun = function(newdata,times) {Lambda1(newdata,times) + Lambda2(newdata,times)}
##     ## Test:
##     B_fun = function(newdata,times) {Lambda1_(newdata,times) + Lambda2_(newdata,times)}
##     out = integral2(A = A_fun,
##                     B = B_fun,
##                     C = C_fun,
##                     D = Lambda1,
##                     end = t,
##                     data = data,
##                     jump_points = jump_points,
##                     chunks = chunks,
##                     integrate = FALSE)
##     ## Alternative form:
##     h_fun = function(newdata,times) exp(-(Lambda1_(newdata,times) + Lambda2_(newdata,times)))
##     gamma_fun = function(newdata,times) {
##         exp(Lambda1_(newdata,times) + Lambda2_(newdata,times) + Gamma_(newdata,times))*
##             at_risk_fun(newdata, times)
##     }
##     bar_fun = function(newdata,times) {Lambda1(newdata,times) + Lambda2(newdata,times)}
##     term_plus1 = integral1(A = h_fun,B = Lambda1,end = t,data = data,jump_points = jump_points,chunks = chunks,integrate = FALSE)    
##     term_plus2 = integral1(A = gamma_fun,B = bar_fun,end = t,data = data,jump_points = jump_points,chunks = chunks,integrate = FALSE)
##     term_minus = integral2(A = h_fun,B = Lambda1,C = gamma_fun,D = bar_fun,end = t,data = data,jump_points = jump_points,chunks = chunks,integrate = FALSE)
##     return(list(out = out, term_plus1 = term_plus1, term_plus2 = term_plus2, term_minus = term_minus))
## }

## at_risk_fun <- function(time_name = "time"){
##     out = function(newdata,times){
##         spring_times = newdata[[time_name]]
##         grid_mat = matrix(times,ncol = length(times),nrow = length(spring_times),byrow = TRUE)
##         event_mat = apply(grid_mat, 2, function(x) 1*(spring_times >= x))
##         return(event_mat)
##     }
##     return(out)
## }



######################################################################
### martingale-int.R ends here
