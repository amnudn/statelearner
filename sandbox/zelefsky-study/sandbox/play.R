### play.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 18 2023 (15:39) 
## Version: 
## Last-Updated: Jul 18 2023 (16:08) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(lava)
library(data.table)
library(prodlim)
m <- lvm()
formula_event <- ~f(age,0)
formula_cens <- ~f(age,0)
event_scale <- 1/1000
cens_scale <- 1/1000
lava::distribution(m, ~age) <- lava::normal.lvm(mean = 60,sd = 11)
lava::distribution(m, ~biomarker) <- lava::normal.lvm(mean = 0,sd = 1)
lava::regression(m) <- eventtime ~ 0.01*age+0.4*sex
lava::regression(m) <- censtime ~ 0.03*age-0.4*sex
lava::regression(m) <- S_event ~ exp(0.4*como-0.5*sex)
lava::regression(m) <- S_cens ~ exp(-.4*como+0.8*sex)
lava::distribution(m,~sex+como) <- lava::binomial.lvm()
lava::distribution(m,~censtime) <- lava::coxWeibull.lvm(scale=event_scale,shape = ~S_cens)
lava::distribution(m,~eventtime) <- lava::coxWeibull.lvm(scale=cens_scale,shape=~S_event)
m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0), "event")
## lava::regression(m) <- stats::update(formula_event, "eventtime~.")
## lava::regression(m) <- stats::update(formula_cens, "censtime~.")
d <- setDT(sim(m,1000))
d[,table(event)]
pdf("~/tmp/test.pdf",width = 12,height = 6)
par(mfrow = c(1,2))
plot(prodlim(Hist(time,event)~sex+age,data = d),xlim = c(0,50),newdata = expand.grid(age = c(20,60,90),sex =c(0,1)),type = "risk",plot.main = "Event probability",confint = FALSE,legend.cex = 0.5)
plot(prodlim(Hist(time,event)~sex+age,data = d,rev = TRUE),xlim = c(0,50),newdata = expand.grid(age = c(20,60,90),sex =c(0,1)),plot.main = "Censoring probability",confint = FALSE,legend.cex = 0.5)
dev.off()
######################################################################
### play.R ends here
