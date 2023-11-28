### sim_zelefsky.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov 12 2023 (15:55) 
## Version: 
## Last-Updated: Nov 28 2023 (09:57) 
##           By: Anders Munch
##     Update #: 74
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(survival)
library(data.table)
library(lava)
library(here)

if(FALSE){
    zelefsky = setDT(get(load(here("data/zelefsky.rda"))))
    zelefsky$status=as.numeric(as.character(factor(zelefsky$recur,levels=c("No","Yes"),labels=c("0","1"))))
    zelefsky$logPSA <- log(zelefsky$psa)
    zelefsky$sDose <- as.vector(scale(zelefsky$dose))
    zelefsky$hormonesYes <- zelefsky$hormones=="Yes"
    for (v in levels(zelefsky$stage)){    
        zelefsky[,paste("stage",v,sep="")] <- as.numeric(zelefsky$stage==v)
    }
    zelefsky[]

    zWeib.surv <- survreg(formula = Surv(dmos, recur != "No") ~ log(psa) + stage + ggtot + sDose + hormones,data = zelefsky)
    # censoring hazards
    zWeib.cens <- survreg(formula = Surv(dmos, recur == "No") ~ log(psa) + stage + ggtot + sDose + hormones, data = zelefsky)
    zelefsky_sim_input <- list(surv = list(coef = coef(zWeib.surv),
                                           scale = zWeib.surv$scale),
                               cens = list(coef = coef(zWeib.cens),
                                           scale = zWeib.cens$scale),
                               zelefsky = zelefsky)

    zelefsky_summary <- zelefsky_sim_input
    zelefsky_summary$zelefsky <- NULL
    zelefsky_summary$logPSA <- list(mean = mean(log(zelefsky$psa)), sd = sd(log(zelefsky$psa)))
    zelefsky_summary$ggtot <- list(mean = mean(zelefsky$ggtot), sd = sd(zelefsky$ggtot))
    zelefsky_summary$sDose <- list(mean = mean(zelefsky$sDose), sd = sd(zelefsky$sDose))
    zelefsky_summary$hormones <- list(p = mean(zelefsky$hormonesYes))
    zelefsky_summary$stage$p <- table(zelefsky$stage)[-6]/nrow(zelefsky)
    saveRDS(zelefsky_summary, file = here("R-code", "zelefsky-sim-summary.rda"))
}

zelefsky_summary <- readRDS(here("R-code", "zelefsky-sim-summary.rda"))

simZelefsky <- function(n,censoring = TRUE,simulation_input = zelefsky_summary, ...){
    m <- lvm()
    input.surv = simulation_input$surv
    input.cens = simulation_input$cens
    distribution(m,~time.recur) <- coxWeibull.lvm(scale=exp(-input.surv$coef["(Intercept)"]/input.surv$scale),shape=1/input.surv$scale)
    distribution(m,~time.cens) <- coxWeibull.lvm(scale=exp(-input.cens$coef["(Intercept)"]/input.cens$scale),shape=(1/input.cens$scale))
    ## put code to generate covariates: psa, stage, ggtot, dose100, hormones
    ## alike the real data distributions, maybe even with
    ## regression statements between them to make them correlated

    distribution(m,~logPSA) <- normal.lvm(mean=simulation_input$logPSA$mean,sd=simulation_input$logPSA$sd)
    m <- categorical(m,~stage, K=6, p=c(simulation_input$stage$p), labels=c("T1c","T2a","T2b","T2c","T3ab","T3c"))
    transform(m,stageT2a~stage) <- function(x){1*(x=="T2a")}
    transform(m,stageT2b~stage) <- function(x){1*(x=="T2b")}
    transform(m,stageT2c~stage) <- function(x){1*(x=="T2c")}
    transform(m,stageT3ab~stage) <- function(x){1*(x=="T3ab")}
    transform(m,stageT3c~stage) <- function(x){1*(x=="T3c")}
    distribution(m,~ggtot) <- normal.lvm(mean=simulation_input$ggtot$mean,sd=simulation_input$ggtot$sd)
    distribution(m,~sDose) <- normal.lvm(mean=simulation_input$sDose$mean,sd=simulation_input$sDose$sd)
    
    distribution(m,~hormones) <- binomial.lvm(p=simulation_input$hormones$p)
    regression(m,time.recur~logPSA + stageT2a + stageT2b + stageT2c + stageT3ab + stageT3c + ggtot + sDose + hormones) <- -input.surv$coef[-1]/input.surv$scale
    regression(m,time.cens~logPSA + stageT2a + stageT2b + stageT2c + stageT3ab + stageT3c + ggtot + sDose + hormones) <- -input.cens$coef[-1]/input.cens$scale
    
    if (censoring){
        m <- lava::eventTime(m, dmos ~ min(time.recur = 1, time.cens=0),eventName = "recur")
        d <- sim(m,n)
    } else {
        d <- sim(m,n)
        d$dmos <- d$time.recur
        d$recur <- 1
    }

    ## Something is strange here?
    for (x in 1:ncol(d)){
        d[,x] <- c(d[,x])
    }
    setDT(d)
    d[,psa:=exp(logPSA)]
    d[,hormones:=factor(hormones, levels=c("0","1"), labels = c("No", "Yes"))]
    d[,status:=1*recur==1]
    d[,recur:=factor(recur, levels=c("0","1"), labels = c("No", "Yes"))]
    d[,]
    d
}

## simZelefsky <- function(n,censoring = TRUE,simulation_input = zelefsky_sim_input, ...){
##     m <- lvm()
##     input.surv = simulation_input$surv
##     input.cens = simulation_input$cens
##     zelefsky = simulation_input$zelefsky
##     distribution(m,~time.recur) <- coxWeibull.lvm(scale=exp(-input.surv$coef["(Intercept)"]/input.surv$scale),shape=1/input.surv$scale)
##     distribution(m,~time.cens) <- coxWeibull.lvm(scale=exp(-input.cens$coef["(Intercept)"]/input.cens$scale),shape=(1/input.cens$scale))
##     ## put code to generate covariates: psa, stage, ggtot, dose100, hormones
##     ## alike the real data distributions, maybe even with
##     ## regression statements between them to make them correlated

##     distribution(m,~logPSA) <- normal.lvm(mean=mean(log(zelefsky$psa)),sd=sd(log(zelefsky$psa)))
##     m <- categorical(m,~stage, K=6, p=c(table(zelefsky$stage)[-6]/nrow(zelefsky)), labels=c("T1c","T2a","T2b","T2c","T3ab","T3c"))
##     transform(m,stageT2a~stage) <- function(x){1*(x=="T2a")}
##     transform(m,stageT2b~stage) <- function(x){1*(x=="T2b")}
##     transform(m,stageT2c~stage) <- function(x){1*(x=="T2c")}
##     transform(m,stageT3ab~stage) <- function(x){1*(x=="T3ab")}
##     transform(m,stageT3c~stage) <- function(x){1*(x=="T3c")}
##     distribution(m,~ggtot) <- normal.lvm(mean=mean(zelefsky$ggtot),sd=sd(zelefsky$ggtot))
##     distribution(m,~sDose) <- normal.lvm(mean=mean(zelefsky$sDose),sd=sd(zelefsky$sDose))
    
##     distribution(m,~hormones) <- binomial.lvm(p=mean(zelefsky$hormonesYes))
##     regression(m,time.recur~logPSA + stageT2a + stageT2b + stageT2c + stageT3ab + stageT3c + ggtot + sDose + hormones) <- -input.surv$coef[-1]/input.surv$scale
##     regression(m,time.cens~logPSA + stageT2a + stageT2b + stageT2c + stageT3ab + stageT3c + ggtot + sDose + hormones) <- -input.cens$coef[-1]/input.cens$scale
    
##     if (censoring){
##         m <- lava::eventTime(m, dmos ~ min(time.recur = 1, time.cens=0),eventName = "recur")
##         d <- sim(m,n)
##     } else {
##         d <- sim(m,n)
##         d$dmos <- d$time.recur
##         d$recur <- 1
##     }

##     ## Something is strange here?
##     for (x in 1:ncol(d)){
##         d[,x] <- c(d[,x])
##     }
##     setDT(d)
##     d[,psa:=exp(logPSA)]
##     d[,hormones:=factor(hormones, levels=c("0","1"), labels = c("No", "Yes"))]
##     d[,status:=1*recur==1]
##     d[,recur:=factor(recur, levels=c("0","1"), labels = c("No", "Yes"))]
##     d[,]
##     d
## }

simZelefsky_wrapper <- function(n,censoring = TRUE,simulation_input = zelefsky_summary, ...){
    raw_dat = simZelefsky(n = n, censoring = censoring, simulation_input = simulation_input, ...)
    out = raw_dat[, .(X1 = logPSA, X2 = stage, X3 = ggtot, X4 = sDose, X5 = hormones, time = dmos, status = 1*status, true_time = time.recur, cens_time = time.cens)]
    return(out)
}

simZelefsky_indep_cens_wrapper <- function(n,simulation_input = zelefsky_summary, ...){
    indep_cens = c("(Intercept)" = 4, rep(0, length(simulation_input$cens$coef)-1))
    ## indep_cens = c(simulation_input$cens$coef["(Intercept)"], rep(0, length(simulation_input$cens$coef)-1))
    simulation_input$cens$coef = indep_cens
    raw_dat = simZelefsky(n = n, censoring = TRUE, simulation_input = simulation_input, ...)
    out = raw_dat[, .(X1 = logPSA, X2 = stage, X3 = ggtot, X4 = sDose, X5 = hormones, time = dmos, status = 1*status, true_time = time.recur, cens_time = time.cens)]
    return(out)
}

simZelefsky_simple_effect_wrapper <- function(n,simulation_input = zelefsky_summary, ...){
    simple_effect = c("(Intercept)" = 5, "log(psa)" = -0.5, rep(0, length(simulation_input$cens$coef)-2))
    ## indep_cens = c(simulation_input$cens$coef["(Intercept)"], rep(0, length(simulation_input$cens$coef)-1))
    simulation_input$surv$coef = simple_effect
    raw_dat = simZelefsky(n = n, censoring = TRUE, simulation_input = simulation_input, ...)
    out = raw_dat[, .(X1 = logPSA, X2 = stage, X3 = ggtot, X4 = sDose, X5 = hormones, time = dmos, status = 1*status, true_time = time.recur, cens_time = time.cens)]
    return(out)
}

simZelefsky_noise_wrapper <- function(n,censoring = TRUE,simulation_input = zelefsky_summary, p = 5,...){
    raw_dat = simZelefsky(n = n, censoring = censoring, simulation_input = simulation_input, ...)
    out = raw_dat[, .(X1 = logPSA, X2 = stage, X3 = ggtot, X4 = sDose, X5 = hormones, time = dmos, status = 1*status, true_time = time.recur, cens_time = time.cens)]
    ## Add additional noise
    if(p>0){
        Xp = as.data.table(matrix(rnorm(n*p), n, p))
        names(Xp) = paste0("X", 5+(1:p))
        out = cbind(Xp, out)
    }
    return(out)
}

## ## Settings from Zelefsky:
## surv$coef <-  c(7.6627958, -0.5111233, -0.6097065, -0.9189468, -1.1148973, -1.0502660, -1.1500966, -0.1302955,  0.3786362,  0.2308239)

## ## surv
## c(7.7,-0.5,-0.6,-0.9,-1.1,-1.1,-1.2,-0.1,0.4,0.2)
## 0.9

## ## cens
## c(3.3,0.1,-0.1,0.2,0.4,0.5,0.5,0.0,-0.5,-0.2)
## 0.6


######################################################################
### sim_zelefsky.R ends here
