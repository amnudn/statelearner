
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
