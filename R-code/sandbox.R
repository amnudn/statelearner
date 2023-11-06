### sandbox.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Oct 27 2023 (15:40) 
## Version: 
## Last-Updated: Nov  2 2023 (12:06) 
##           By: Anders Munch
##     Update #: 216
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(targets)
library(riskRegression)
library(survival)
library(prodlim)
tar_source("./functions/")

## Sandbox
## Get something to play with


## Reproduce CSC functionality first:
set.seed(44)
d <- sampleData(300,outcome="competing.risks")
aa <- CSC(list(Hist(time,event)~X8, Hist(time,event)~X1+X2+X8), data=d)
d[, state1 := 1*(event == 1)]
d[, state2 := 1*(event == 2)]
cox1 <- coxph(Surv(time,event = state1)~X8, data=d, x = TRUE, y = TRUE)
cox2 <- coxph(Surv(time,state2)~X8, data=d, x = TRUE, y = TRUE)
time_jumps <- d[time <= .5, sort(time)]
## time_jumps <- c(0,d[time <= .5, sort(time)])
cumHazard1 <- predictCox(cox1, newdata = d[1:5], times = time_jumps, type = "cumhazard")$cumhazard
cumHazard2 <- predictCox(cox2, newdata = d[1:5], times = time_jumps, type = "cumhazard")$cumhazard
Hazard1 <- predictCox(cox1, newdata = d[1:5], times = time_jumps, type = "hazard")$hazard
Hazard2 <- predictCox(cox2, newdata = d[1:5], times = time_jumps, type = "hazard")$hazard
## Slight differences are explained by difference between S(t) and S(t-)
as.numeric(predictRisk(aa, newdata = d[1:5], times = .5, cause = 1, product.limit = 0))
## rowSums(exp(- cumHazard1-cumHazard2)*Hazard1)
rowSums(cbind(1,exp(- cumHazard1[,-ncol(cumHazard1), drop = FALSE]-cumHazard2[,-ncol(cumHazard2), drop = FALSE]))*Hazard1)

## Testing:
## Example functions:
fun1 <- function(newdata, times){
    predictCox(cox1, newdata = newdata, times = times, type = "cumhazard")$cumhazard
}
## Example of function using t-:
min_time_diff <- d[, min(diff(sort(time)))]
fun2 <- function(newdata, times){
    times0 = times-(min_time_diff/2)
    times0 = pmax(0, times0)
    ch1 = predictCox(cox1, newdata = newdata, times = times0, type = "cumhazard")$cumhazard
    ## return(ch1)
    ch2 = predictCox(cox2, newdata = newdata, times = times0, type = "cumhazard")$cumhazard
    return(-ch1-ch2)
}
## These are the same! -- but last expression does not agree with predict...?
internal_integral1(A = fun2, C = fun1, end = .5, jump_points = d[, time], newdata = d[1:5], integrate = TRUE)
rowSums(cbind(1,exp(- cumHazard1[,-ncol(cumHazard1), drop = FALSE]-cumHazard2[,-ncol(cumHazard2), drop = FALSE]))*Hazard1)
## But does not agree with predict?
## as.numeric(predictRisk(aa, newdata = d[1:5], times = .5, cause = 1, product.limit = 0))

## Consistency:
internal_integral1(A = fun2, C = fun1, end = .5, jump_points = d[, time], newdata = d[1:5], integrate = TRUE)
internal_integral1(A = fun2, C = fun1, end = .25, jump_points = d[, time], newdata = d[1:5], integrate = TRUE) +
    internal_integral1(A = fun2, C = fun1, start = .25, end = .5, jump_points = d[, time], newdata = d[1:5], integrate = TRUE)
internal_integral1(A = fun2, C = fun1, end = d[time <= .5 & event == 1, sort(time)[3]], jump_points = d[, time], newdata = d[1:5], integrate = TRUE) +
    internal_integral1(A = fun2, C = fun1, start = d[time <= .5 & event == 1, sort(time)[3]], end = .5, jump_points = d[, time], newdata = d[1:5], integrate = TRUE)
internal_integral1(A = fun2, C = fun1, end = d[time <= .5 & event == 1, sort(time)[3]]-min_time_diff, jump_points = d[, time], newdata = d[1:5], integrate = TRUE) +
    internal_integral1(A = fun2, C = fun1, start = d[time <= .5 & event == 1, sort(time)[3]]-min_time_diff, end = .5, jump_points = d[, time], newdata = d[1:5], integrate = TRUE)


internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = .5, jump_points = d[, time], data = d[1:5], integrate = TRUE)
internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = .25, jump_points = d[, time], data = d[1:5], integrate = TRUE) +
    internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, start = .25, end = .5, jump_points = d[, time], data = d[1:5], integrate = TRUE)
internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = d[time <= .5 & event == 1, sort(time)[3]], jump_points = d[, time], data = d[1:5], integrate = TRUE) +
    internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, start = d[time <= .5 & event == 1, sort(time)[3]], end = .5, jump_points = d[, time], data = d[1:5], integrate = TRUE)
internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = d[time <= .5 & event == 1, sort(time)[3]]-min_time_diff, jump_points = d[, time], data = d[1:5], integrate = TRUE) +
    internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, start = d[time <= .5 & event == 1, sort(time)[3]]-min_time_diff, end = .5, jump_points = d[, time], data = d[1:5], integrate = TRUE)


internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = d[, sort(time)[10]], jump_points = d[, time], data = d[1:5], integrate = TRUE)
internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = d[, sort(time)[9]], jump_points = d[, time], data = d[1:5], integrate = TRUE) +
    internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, start = d[, sort(time)[9]], end = d[, sort(time)[10]], jump_points = d[, time], data = d[1:5], integrate = TRUE)

system.time(internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE))
system.time(integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE))
system.time(integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE, chunks = 1))

internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d[1:5], integrate = TRUE)
integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d[1:5], integrate = TRUE)
integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d[1:5], integrate = TRUE, chunks = 5)

library(rbenchmark)

benchmark("internal" =  internal_integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE),
          "wrapper" =  integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE),
          "wrapper_chunk2" =  integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE, chunks = 2),
          "wrapper_chunk5" =  integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE, chunks = 5),
          "wrapper_chunk10" =  integral1(A = function(newdata,times) exp(fun2(newdata,times)), B = fun1, end = 2, jump_points = d[, time], data = d, integrate = TRUE, chunks = 10),
          replications = 200)


99

## Function on specific form:
## Input: (newdata, times)
## Output: matrix of CH? (row for each person in new data, column for each time point)

## Testing consistency of nested integral function:
## A=1,B=fun2,C=fun1,D=fun1
## Analytically, we have that
## \int_start^end A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
## = \int_start^m A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
##   + \int_m^end A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
## = \int_start^m A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
##   + \int_m^end A(t) [\int_0^m exp(B(s)) C(ds)] D(dt)
##   + \int_m^end A(t) [\int_m^t exp(B(s)) C(ds)] D(dt)
## = \int_start^m A(t) [\int_0^t exp(B(s)) C(ds)] D(dt)
##   + [\int_m^end A(t) D(dt)] * [\int_0^m exp(B(s)) C(ds)] 
##   + \int_m^end A(t) [\int_m^t exp(B(s)) C(ds)] D(dt)
## Confirm:

internal_integral2(A = function(newdata,times) exp(fun2(newdata,times)),
                   B = fun1,
                   C = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)),
                   D = fun1,
                   end = 2,
                   jump_points = d[, time],
                   data = d[1:5])
internal_integral2(A = function(newdata,times) exp(fun2(newdata,times)),
                   B = fun1,
                   C = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)),
                   D = fun1,
                   end = 1,jump_points = d[, time],data = d[1:5]) +
    integral1(A = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)),
              B = fun1,
              start = 1, end = 2,
              jump_points = d[, time],data = d[1:5]) *
    integral1(A = function(newdata,times) exp(fun2(newdata,times)),
              B = fun1,
              start = 0, end = 1,
              jump_points = d[, time],data = d[1:5]) +
    internal_integral2(A = function(newdata,times) exp(fun2(newdata,times)),
                       B = fun1,
                       C = function(newdata,times) matrix(1,nrow = nrow(newdata),ncol = length(times)),
                       D = fun1,start = 1,end = 2,jump_points = d[, time],data = d[1:5])


internal_integral2(A = function(newdata, times) matrix(1,nrow = nrow(newdata),ncol = length(times)),                   
                   B = function(newdata, times) matrix(c(rep(1, nrow(newdata)), rep(2, (length(times)-1)*nrow(newdata))),nrow = nrow(newdata),ncol = length(times)),
                   C = function(newdata,times) exp(-fun1(newdata,times)),
                   D = fun1,
                   end = 2,
                   jump_points = d[, time],
                   data = d[1:5])
internal_integral1(A = function(newdata,times) exp(-fun1(newdata,times)),
                   B = fun1,
                   end = 2,
                   jump_points = d[, time],
                   data = d[1:5])

internal_integral2(A = function(newdata,times) exp(fun2(newdata,times)),
                   B = fun1,
                   C = function(newdata,times) exp(-fun1(newdata,times)),
                   D = fun1,
                   end = 2,
                   jump_points = d[, time],
                   data = d[1:5])

integral2(A = function(newdata,times) exp(fun2(newdata,times)),
          B = fun1,
          C = function(newdata,times) exp(-fun1(newdata,times)),
          D = fun1,
          end = 2,
          jump_points = d[, time],
          data = d[1:5])
integral2(A = function(newdata,times) exp(fun2(newdata,times)),
          B = fun1,
          C = function(newdata,times) exp(-fun1(newdata,times)),
          D = fun1,
          end = 1,jump_points = d[, time],data = d[1:5]) +
    integral1(A = function(newdata,times) exp(-fun1(newdata,times)),
              B = fun1,
              start = 1, end = 2,
              jump_points = d[, time],data = d[1:5]) *
    integral1(A = function(newdata,times) exp(fun2(newdata,times)),
              B = fun1,
              start = 0, end = 1,
              jump_points = d[, time],data = d[1:5]) +
    integral2(A = function(newdata,times) exp(fun2(newdata,times)),
              B = fun1,
              C = function(newdata,times) exp(-fun1(newdata,times)),
              D = fun1,
              start = 1,end = 2,jump_points = d[, time],data = d[1:5])


benchmark("internal" = internal_integral2(A = function(newdata,times) exp(fun2(newdata,times)),B = fun1,C = function(newdata,times) exp(-fun1(newdata,times)),D = fun1,end = 8,jump_points = d[, time],data = d),
          "wrapper" = integral2(A = function(newdata,times) exp(fun2(newdata,times)),B = fun1,C = function(newdata,times) exp(-fun1(newdata,times)),D = fun1,end = 8,jump_points = d[, time],data = d),
          "wrapper_chunk2" = integral2(A = function(newdata,times) exp(fun2(newdata,times)),B = fun1,C = function(newdata,times) exp(-fun1(newdata,times)),D = fun1,end = 8,jump_points = d[, time],data = d, chunks = 2),
          "wrapper_chunk5" = integral2(A = function(newdata,times) exp(fun2(newdata,times)),B = fun1,C = function(newdata,times) exp(-fun1(newdata,times)),D = fun1,end = 8,jump_points = d[, time],data = d, chunks = 5),
          "wrapper_chunk10" = integral2(A = function(newdata,times) exp(fun2(newdata,times)),B = fun1,C = function(newdata,times) exp(-fun1(newdata,times)),D = fun1,end = 8,jump_points = d[, time],data = d, chunks = 10),
          replications = 200)


## Simulate some data
source("./functions/hely-colon-sim.R")
#-- fit models to data; 
fit.colon.cr <- fit.colon.fun(      
    formula.1=Surv(time, event==1)~rx+sex+differ+age+nodes.squared+obstruct+perfor+adhere+extent+surg+rx*sex+rx*perfor+rx*age,
    formula.2=Surv(time, event==2)~rx+sex+nodes+differ+age+obstruct+adhere+extent+surg,
    formula.0=Surv(time, event==0)~rx+sex+nodes+differ+age+obstruct+perfor+adhere+extent+surg,
    formula.treat=rx~sex+age+nodes+differ+obstruct+perfor+adhere+extent+surg,
    d=colon.cr)
#-- here we simulate data (fixed seed); 
set.seed(31)      
sim.colon.cr <- synthesize.colon.fun(fit.colon=fit.colon.cr,
                                     d=colon.cr,
                                     name.treat="rx",
                                     event.name="event")
## Setup data on right format
setnames(sim.colon.cr, "event", "status")
sim.colon.cr[, A := 1*(rx.num>1)]

## setup some models:
cause1_cox <- coxph(Surv(time, status==1)~A+sex+differ+age+obstruct+perfor+adhere,
                    data = sim.colon.cr,
                    x = TRUE, y = TRUE)
cause2_cox <- coxph(Surv(time, status==2)~A+sex,
                    data = sim.colon.cr,
                    x = TRUE, y = TRUE)
cens_cox <- coxph(Surv(time, status==0)~A+sex+nodes+differ+age+obstruct+perfor,
                  data = sim.colon.cr,
                  x = TRUE, y = TRUE)

treat_model <- glm(A ~ sex+nodes+differ+age+obstruct+perfor, family = binomial(), data = sim.colon.cr)

L1_test <- construct_pred_fun(cause1_cox)
L2_test <- construct_pred_fun(cause2_cox)
L0_test <- construct_pred_fun(cens_cox)
pi_test <- construct_pred_fun(treat_model)

## >1 x >1
L0_test(sim.colon.cr[1:5], times = seq(100, 500, 100))
pi_test(sim.colon.cr[1:5])
at_risk_fun(sim.colon.cr[1:5], times = seq(200, 500, 100))
## > 1 x =1
L0_test(sim.colon.cr[1:5], times = 100)
at_risk_fun(sim.colon.cr[1:5], times = 100) 
## =1 x >1
L0_test(sim.colon.cr[1], times = seq(100, 500, 100))
pi_test(sim.colon.cr[1])
at_risk_fun(sim.colon.cr[1], times = seq(100, 500, 100)) 
## =1 x =1
L0_test(sim.colon.cr[1], times = 100)
at_risk_fun(sim.colon.cr[1], times = 100) 

L0_test(sim.colon.cr[1:5], times = 100)






internal_integral1(A = function(newdata,times) matrix(0,nrow = nrow(newdata),ncol = length(times)),
                    B = at_risk_all,
                    C = L1_test,
                    newdata = sim.colon.cr[1:4],
                    end = 2700,
                    jump_points = sort(unique(sim.colon.cr[, time])))

internal_integral1(A = L2_test,C = L1_test,newdata = sim.colon.cr[1:4],end = 2700,jump_points = sort(unique(sim.colon.cr[, time])))

internal_integral1(A = L2_test,B = at_risk_all,C = L1_test,newdata = sim.colon.cr[1:4],end = 2700,jump_points = sort(unique(sim.colon.cr[, time]))) +
    internal_integral1(A = L2_test,B = function(newdata, times) {1 - at_risk_all(newdata,times)},C = L1_test,newdata = sim.colon.cr[1:4],end = 2700,jump_points = sort(unique(sim.colon.cr[, time])))

diag(predictCox(cause1_cox, newdata = sim.colon.cr[1:4], times = sim.colon.cr[1:4, time], type = "cumhazard")$cumhazard)

diag(L1_test(newdata = sim.colon.cr[1:4], times = sim.colon.cr[1:4, time]))


naiv_est(data = sim.colon.cr,
          t = 200,
          Lambda1 = L1_test,
          Lambda2 = L2_test,
          jump_points = sort(unique(sim.colon.cr[, time])))

naiv_est(data = sim.colon.cr,
          t = 200,
          Lambda1 = L1_test,
          Lambda2 = L2_test,
          jump_points = sort(unique(sim.colon.cr[, time])),
          chunks = 5)



f1_test(sim.colon.cr[c(1,2,3:7)], c(100, 200))

naiv_est(data = sim.colon.cr[c(2,1,3:7)],
         t = 200,
         Lambda1 = L1_test,
         Lambda2 = L2_test,
         jump_points = sort(unique(sim.colon.cr[, time])),
         collapse = FALSE)

Q1_fun(data = sim.colon.cr,
       t = 200,
       Lambda1 = L1_test,
       Lambda2 = L2_test,
       Gamma = L0_test,
       jump_points = sort(unique(sim.colon.cr[, time])),
       collapse = TRUE,
       check = FALSE)

Q1_fun(data = sim.colon.cr,
       t = 200,
       Lambda1 = L1_test,
       Lambda2 = L2_test,
       Gamma = L0_test,
       jump_points = sort(unique(sim.colon.cr[, time])),
       collapse = TRUE) -
    Q2_fun(data = sim.colon.cr,
           t = 200,
           Lambda1 = L1_test,
           Lambda2 = L2_test,
           Gamma = L0_test,
           jump_points = sort(unique(sim.colon.cr[, time])),
           collapse = TRUE)



R1_fun(data = sim.colon.cr,
       t = 200,
       Lambda1 = L1_test,
       Lambda2 = L2_test,
       Gamma = L0_test,
       jump_points = sort(unique(sim.colon.cr[, time])),
       collapse = TRUE) -
    R2_fun(data = sim.colon.cr,
           t = 200,
           Lambda1 = L1_test,
           Lambda2 = L2_test,
           Gamma = L0_test,
           jump_points = sort(unique(sim.colon.cr[, time])),
           collapse = TRUE) 

R1_fun(data = sim.colon.cr,
           t = 200,
           Lambda1 = L1_test,
           Lambda2 = L2_test,
           Gamma = L0_test,
           jump_points = sort(unique(sim.colon.cr[, time])),
       collapse = TRUE)

R1_count_fun(data = sim.colon.cr,
             t = 200,
             Lambda1 = L1_test,
             Lambda2 = L2_test,
             Gamma = L0_test,
             jump_points = sort(unique(sim.colon.cr[, time])),
             collapse = TRUE) -
    R1_comp_fun(data = sim.colon.cr,
                t = 200,
                Lambda1 = L1_test,
                Lambda2 = L2_test,
                Gamma = L0_test,
                jump_points = sort(unique(sim.colon.cr[, time])),
                collapse = TRUE)

## This is OK!
## R1_count_fun(data = sim.colon.cr,
##              t = 300,
##              Lambda1 = L1_test,
##              Lambda2 = L2_test,
##              Gamma = L0_test,
##              jump_points = sort(unique(sim.colon.cr[, time])),
##              collapse = TRUE) -
##     R1_comp_fun(data = sim.colon.cr,
##                 t = 300,
##                 Lambda1 = L1_test,
##                 Lambda2 = L2_test,
##                 Gamma = L0_test,
##                 jump_points = sort(unique(sim.colon.cr[, time])),
##                 collapse = TRUE) -
##     (R2_count_fun(data = sim.colon.cr,
##                   t = 300,
##                   Lambda1 = L1_test,
##                   Lambda2 = L2_test,
##                   Gamma = L0_test,
##                   jump_points = sort(unique(sim.colon.cr[, time])),
##                   collapse = TRUE) -
##      R2_comp_fun(data = sim.colon.cr,
##                  t = 300,
##                  Lambda1 = L1_test,
##                  Lambda2 = L2_test,
##                  Gamma = L0_test,
##                  jump_points = sort(unique(sim.colon.cr[, time])),                
##                  collapse = TRUE))
## Q1_fun(data = sim.colon.cr,
##        t = 300,
##        Lambda1 = L1_test,
##        Lambda2 = L2_test,
##        Gamma = L0_test,
##        jump_points = sort(unique(sim.colon.cr[, time])),
##        collapse = TRUE) -
##     Q2_fun(data = sim.colon.cr,
##            t = 300,
##            Lambda1 = L1_test,
##            Lambda2 = L2_test,
##            Gamma = L0_test,
##            jump_points = sort(unique(sim.colon.cr[, time])),
##            collapse = TRUE)

mean(termB(data = sim.colon.cr,
           t = 300,
           Lambda1 = L1_test,
           Lambda2 = L2_test,
           Gamma = L0_test,
           jump_points = sort(unique(sim.colon.cr[, time])))-
     termC(data = sim.colon.cr,
           t = 300,
           Lambda1 = L1_test,
           Lambda2 = L2_test,
           Gamma = L0_test,
           jump_points = sort(unique(sim.colon.cr[, time]))))

tt <- one_step(data = sim.colon.cr,
               t = 300,
               Lambda1 = L1_test,
               Lambda2 = L2_test,
               Gamma = L0_test,
               pi = pi_test,
               jump_points = sort(unique(sim.colon.cr[, time])), collapse = FALSE)

mean(tt$naiv_i)
mean(tt$W_i*tt$A_i)
mean(tt$W_i*tt$B_i)
mean(tt$W_i*tt$C_i)

Q2_fun(data = sim.colon.cr,
       t = 300,
       Lambda1 = L1_test,
       Lambda2 = L2_test,
       Gamma = L0_test,
       jump_points = sort(unique(sim.colon.cr[, time])),
       collapse = FALSE)[1:5]
## Q2_fun_v2(data = sim.colon.cr,
##        t = 300,
##        Lambda1 = L1_test,
##        Lambda2 = L2_test,
##        Gamma = L0_test,
##        jump_points = sort(unique(sim.colon.cr[, time])),
##        collapse = FALSE)[1:5]

tt <- Q2_check(data = sim.colon.cr[1:5],
               t = 300,
               Lambda1 = L1_test,
               Lambda2 = L2_test,
               Gamma = L0_test,
               jump_points = sort(unique(sim.colon.cr[, time])))

tt$out[, 1:6]
tt$term_plus1[, 1:6]
tt$term_plus2[, 1:6]
tt$term_minus[, 1:6]

rowSums(tt$out)
rowSums(tt$term_plus1)*rowSums(tt$term_plus2)-rowSums(tt$term_minus)

rowSums(tt$out)
rowSums(tt$term_plus1)*rowSums(tt$term_plus2)-rowSums(cbind(0,tt$term_minus[,-ncol(tt$term_minus)]))

rowSums(tt$out[,1:2])
rowSums(tt$term_plus1[,1:2])*rowSums(tt$term_plus2[,1:2])-rowSums(tt$term_minus[,1:2])

tt$out[,1]
tt$term_plus1[,1]*tt$term_plus2[,1]-tt$term_minus[,1]


tt$out[,1:2]
tt$term_plus1[,1:2]
tt$term_plus2[,1:2]
tt$term_minus[,1:2]



R2_fun(data = sim.colon.cr,
       t = 200,
       Lambda1 = L1_test,
       Lambda2 = L2_test,
       Gamma = L0_test,
       jump_points = sort(unique(sim.colon.cr[, time])),
       collapse = TRUE)
R2_count_fun(data = sim.colon.cr,
             t = 200,
             Lambda1 = L1_test,
             Lambda2 = L2_test,
             Gamma = L0_test,
             jump_points = sort(unique(sim.colon.cr[, time])),
             collapse = TRUE) -
    R2_comp_fun(data = sim.colon.cr,
                t = 200,
                Lambda1 = L1_test,
                Lambda2 = L2_test,
                Gamma = L0_test,
                jump_points = sort(unique(sim.colon.cr[, time])),
                collapse = TRUE)


## OK!

diag(L1_test(sim.colon.cr[1:50], sim.colon.cr[1:50, time]))

integral0_itime(L1_test,times = sim.colon.cr[1:50, time],data = sim.colon.cr[1:50])

integral0_itime(L1_test,times = sim.colon.cr[1:110, time],data = sim.colon.cr[1:110], chunks = 1)


######################################################################
### sandbox.R ends here
