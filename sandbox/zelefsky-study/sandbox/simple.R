library(riskRegression)
library(survival)
library(lava)
library(data.table)
library(prodlim)
m <- lvm()
formula_event <- ~f(age,-0.01)+f(sex,0.3)+f(quad.age,0.01)
formula_cens <- ~f(age,0)+f(sex,1.8)
event_scale <- 1/1000
cens_scale <- 1/1000
lava::distribution(m, ~age) <- lava::normal.lvm(mean = 60,sd = 11)
transform(m,quad.age~age) <- function(x)(x$age-60)^2
lava::distribution(m,~sex) <- lava::binomial.lvm(p = .3)
lava::distribution(m,~censtime) <- lava::coxWeibull.lvm(scale=event_scale)
lava::distribution(m,~eventtime) <- lava::coxWeibull.lvm(scale=cens_scale)
m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0), "event")
lava::regression(m) <- stats::update(formula_event, "eventtime~.")
lava::regression(m) <- stats::update(formula_cens, "censtime~.")
d <- setDT(sim(m,1000))
d[,table(event)]
dlarge <- setDT(sim(m,10000))
dlarge[,dummy := 1]
pdf("~/tmp/test.pdf",width = 12,height = 6)
par(mfrow = c(1,2))
plot(prodlim(Hist(time,event)~sex+age,data = d),xlim = c(0,50),newdata = expand.grid(sex =c(0,1),age = c(20,50,90)),type = "risk",plot.main = "Event probability",confint = FALSE,legend.cex = 0.5)
plot(prodlim(Hist(censtime,dummy)~sex,data = dlarge),xlim = c(0,50),newdata = expand.grid(sex =c(0,1)),plot.main = "Censoring probability",confint = FALSE,legend.cex = 0.5,type = "risk")
plot(prodlim(Hist(time,event)~1,data = d,rev = TRUE),xlim = c(0,50),plot.main = "Censoring probability",confint = FALSE,legend.cex = 0.5,add = TRUE,type = "risk",lty = 3)
dev.off()
d[,sex_age := 1]
d[sex == 1,sex_age := age]
d[,age_sex := 1]
d[sex == 0,age_sex := age]
dlarge[,sex_age := 1]
dlarge[sex == 1,sex_age := age]
dlarge[,age_sex := 1]
dlarge[sex == 0,age_sex := age]
mod1 <- coxph(Surv(time,event)~sex+sex_age+quad.age,data = d,x = TRUE,y = TRUE)
mod2 <- coxph(Surv(time,event)~sex+age_sex+quad.age,data = d,x = TRUE,y = TRUE)
dgm <- coxph(Surv(time,event)~sex+age+quad.age,data = d,x = TRUE,y = TRUE)
x <- Score(list(mod1 = mod1,mod2 = mod2,dgm = dgm),data = dlarge,formula = Hist(time,event)~1,metrics = "brier",se.fit = FALSE,times = 15)
x1 <- Score(list(mod1 = mod1,mod2 = mod2,dgm = dgm),data = dlarge,formula = Hist(time,event)~sex,metrics = "brier",se.fit = FALSE,times = 15)
## x
## x1
setkey(x$Brier$score,Brier)
setkey(x1$Brier$score,Brier)
x$Brier$score
x1$Brier$score
