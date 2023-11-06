### zelefsky-empirical-simulate-data.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov 9 2018 (9:49)
## Version:
#----------------------------------------------------------------------
##
### Commentary: Function simulation data similar to the zelefsky data
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
simZelefsky <- function(n,censoring = TRUE,simulation_input, ...){
    m <- lvm()
    input.surv = simulation_input$surv
    input.cens = simulation_input$cens
    zelefsky = simulation_input$zelefsky
    distribution(m,~time.recur) <- coxWeibull.lvm(scale=exp(-input.surv$coef["(Intercept)"]/input.surv$scale),shape=1/input.surv$scale)
    distribution(m,~time.cens) <- coxWeibull.lvm(scale=exp(-input.cens$coef["(Intercept)"]/input.cens$scale),shape=(1/input.cens$scale))
    ## put code to generate covariates: psa, stage, ggtot, dose100, hormones
    ## alike the real data distributions, maybe even with
    ## regression statements between them to make them correlated

    distribution(m,~logPSA) <- normal.lvm(mean=mean(log(zelefsky$psa)),sd=sd(log(zelefsky$psa)))
    m <- categorical(m,~stage, K=6, p=c(table(zelefsky$stage)[-6]/nrow(zelefsky)), labels=c("T1c","T2a","T2b","T2c","T3ab","T3c"))
    transform(m,stageT2a~stage) <- function(x){1*(x=="T2a")}
    transform(m,stageT2b~stage) <- function(x){1*(x=="T2b")}
    transform(m,stageT2c~stage) <- function(x){1*(x=="T2c")}
    transform(m,stageT3ab~stage) <- function(x){1*(x=="T3ab")}
    transform(m,stageT3c~stage) <- function(x){1*(x=="T3c")}
    distribution(m,~ggtot) <- normal.lvm(mean=mean(zelefsky$ggtot),sd=sd(zelefsky$ggtot))
    distribution(m,~sDose) <- normal.lvm(mean=mean(zelefsky$sDose),sd=sd(zelefsky$sDose))    
#    m <- categorical(m,~gtot, K=9, p=c(table(zelefsky$ggtot)[-9]/nrow(zelefsky)))
#    transform(m,ggtot~gtot) <- function(x) x+2
#    distribution(m,~dose1) <- normal.lvm(mean=7020/100, sd=82/100)
#    distribution(m,~dose2) <- normal.lvm(mean=7560/100, sd=82/100)
#    distribution(m,~dose3) <- normal.lvm(mean=8100/100, sd=82/100)
#    w <- rbinom(n, 2, 0.5)+1
#    transform(m,dose100~dose1+dose2+dose3) <- function(x) ((w==1)*x[1])+((w==2)*x[2])+((w==3)*x[3])
   
    distribution(m,~hormones) <- binomial.lvm(p=mean(zelefsky$hormonesYes))
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

######################################################################
### zelefsky-empirical-simulate-data.R ends here
