library(riskRegression)
library(targets)
library(survival)
library(lava)
library(data.table)
library(prodlim)
try(setwd("~/research/SuperVision/Anders/survival-loss/statelearner/empirical-study/"))
tar_source("functions")
m <- lvm()
# B discriminates only when A == 1
formula_event <- ~f(A,0)+f(B,0)+f(AB,-1.5)+f(AC,-1.3)
formula_cens <- ~f(A,0)+f(B,2)+f(C,0)
event_scale <- 1/1000
cens_scale <- 1/1000
lava::distribution(m,~A) <- lava::binomial.lvm(p = .4)
lava::distribution(m,~B) <- lava::normal.lvm()
lava::distribution(m,~C) <- lava::normal.lvm()
transform(m,AB~A+B) <- prod
transform(m,AC~A+C) <- prod
lava::distribution(m,~censtime) <- lava::coxWeibull.lvm(scale=event_scale)
lava::distribution(m,~eventtime) <- lava::coxWeibull.lvm(scale=cens_scale)
m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0), "event")
lava::regression(m) <- stats::update(formula_event, "eventtime~.")
lava::regression(m) <- stats::update(formula_cens, "censtime~.")
d <- setDT(sim(m,5000))
d[,table(event)]
dlarge <- setDT(sim(m,5000))
dlarge[,dummy := 1]
## pdf("~/tmp/test.pdf",width = 12,height = 6)
## pdf("~/tmp/test.pdf",width = 10,height = 10)
## par(mfrow = c(2,2))
## plot(prodlim(Hist(time,event)~A+B,data = d),xlim = c(0,20),newdata = expand.grid(A =c(0,1),B =c(0,1)),type = "risk",plot.main = "Event probability",confint = FALSE,legend.cex = 0.5)
## plot(prodlim(Hist(time,event)~A+C,data = d),xlim = c(0,20),newdata = expand.grid(A = c(0,1),C = c(0,1)),type = "risk",plot.main = "Event probability",confint = FALSE,legend.cex = 0.5)
## plot(prodlim(Hist(censtime,dummy)~B,data = dlarge),xlim = c(0,20),newdata = expand.grid(B =c(0,1)),plot.main = "Censoring probability",confint = FALSE,legend.cex = 0.5,type = "risk")
## plot(prodlim(Hist(time,event)~1,data = d,rev = TRUE),xlim = c(0,20),plot.main = "Censoring probability",confint = FALSE,legend.cex = 0.5,add = TRUE,type = "risk",lty = 3)
## dev.off()
mod_AB <- coxph(Surv(time,event)~A*B,data = d,x = TRUE,y = TRUE)
mod_AC <- coxph(Surv(time,event)~A*C,data = d,x = TRUE,y = TRUE)
## cmod_AB <- coxph(Surv(time,event==0)~B,data = d,x = TRUE,y = TRUE)
## cmod_AC <- coxph(Surv(time,event==0)~A,data = d,x = TRUE,y = TRUE)
dgm <- coxph(Surv(time,event)~A*(B+C),data = d,x = TRUE,y = TRUE)
ttt <- seq(1,15,1)
x <- Score(list(mod_AB = mod_AB,mod_AC = mod_AC,dgm = dgm),data = dlarge,formula = Hist(time,event)~1,metrics = "brier",se.fit = FALSE,times = ttt,contrasts = FALSE,null.model = TRUE,summary = "ibs")
x1 <- Score(list(mod_AB = mod_AB,mod_AC = mod_AC,dgm = dgm),data = dlarge,formula = Hist(time,event)~A+B,metrics = "brier",se.fit = FALSE,times = ttt,contrasts = FALSE,null.model = TRUE,summary = "ibs")
x2 <- statelearner(learners = list(state = list(mod_AB = mod_AB,mod_AC = mod_AC),censoring = list(cmod_AB = cmod_AB,cmod_AC = cmod_AC)),data = dlarge,B = 1,times = ttt,integrate = FALSE)
setkey(x$Brier$score,IBS)
setkey(x1$Brier$score,IBS)
x$Brier$score[times ==  max(ttt), .(model, Brier = round(100*IBS,2))]
x1$Brier$score[times ==  max(ttt), .(model, Brier = round(100*IBS,2))]
x2
