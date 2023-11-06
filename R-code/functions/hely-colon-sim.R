### hely-colon-sim.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Oct 30 2023 (09:07) 
## Version: 
## Last-Updated: Nov  3 2023 (11:25) 
##           By: Anders Munch
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
library(survival)


fit.categorical <- function(formula,data,pseudo.out=FALSE,...){
    pseudo.data <- get.pseudodata(formula,data)
    # introduce all interactions
    pseudo.formula <- update(formula,"pseudo.outcome~(.)*pseudo.factor")
    # fit pseudo model
    fit <- glm(pseudo.formula, data=pseudo.data[pseudo.data$at.risk==1,], family=binomial())
    if (pseudo.out) {
        pseudo.data$pred <- predict(fit, type="response", newdata=pseudo.data)
        predicted <- cbind(0, do.call("cbind", lapply(unique(pseudo.data$pseudo.aa), function(aa) {
            pseudo.data[pseudo.data$pseudo.aa==aa, "pred"]
        })))
        predicted2 <- do.call("cbind", lapply(2:ncol(predicted), function(jj) {
            apply(1-predicted[, 1:(jj-1), drop=FALSE], 1, prod)*predicted[, jj]
        }))
        colnames(predicted2) <- unique(pseudo.data$pseudo.aa)
    }
    # add the call
    fit$catfit.call <- match.call()
    # add the outcome values
    fit$y <- attr(pseudo.data,"outcome.values")
    # add a new S3 class
    class(fit) <- c("catfit",class(fit))
    # return the fit
    if (pseudo.out) {
        return(list(fit, predicted2))
    } else {
        return(list(fit))
    }
}

get.pseudodata <- function(formula,data,atrisk=TRUE){
    d <- model.frame(formula,data)
    # should check that outcome is categorical 
    N <- NROW(d)
    d$id <- 1:N
    Y <- d[[1]]
    outcome <- names(d)[[1]]
    # HOV: need to deal with factors and factor levels
    # that cannot be coerced to numeric
    if (is.factor(Y))
        outcome.values <- levels(Y)
    else
        outcome.values <- sort(unique(Y))
    pseudo.data <- merge(d,data.frame(id=rep(1:N,each=length(outcome.values)),
                                      pseudo.aa=rep(outcome.values,N)),by="id")
    # should check that no variable in the original data has such name 
    pseudo.data$pseudo.outcome <- as.numeric(pseudo.data[[outcome]]==pseudo.data$pseudo.aa)
    pseudo.data$pseudo.factor <- factor(pseudo.data$pseudo.aa)
    if (atrisk[[1]]==TRUE){
        # indicator of 'at-risk'
        pseudo.data$at.risk <- 0
        pseudo.data$at.risk[pseudo.data[[outcome]]>=pseudo.data$pseudo.aa] <- 1
    }
    # add the outcome name and the unique values
    attr(pseudo.data,"outcome") <- outcome
    attr(pseudo.data,"outcome.values") <- outcome.values
    # return data
    return(pseudo.data)
}

#-- create factor variable out of relevant covariates:
colon$differ <- factor(colon$differ)        
colon$extent <- factor(colon$extent)    

#-- create a dataset with one row per subject:
colon.dt <- merge(setDT(colon)[etype==2], setDT(colon)[etype==1, c("id", "time", "etype", "status"), with=FALSE], by="id")
colon.dt[, time.death:=time.x]  
colon.dt[, time.recurrence:=time.y]     
colon.dt[, status.death:=status.x]       
colon.dt[, status.recurrence:=status.y]
colon.dt[status.death==1 & status.recurrence==0, event:=2]
colon.dt[time.recurrence<=time.death & status.recurrence==1, event:=1]
colon.dt[status.death==0 & status.recurrence==0, event:=0] 
colon.dt <- colon.dt[, !(names(colon.dt) %in% c("status", grep("\\.x|\\.y", names(colon.dt), value=TRUE))), with=FALSE]

#-- create survival data set with event=1 if dead (with/without recurrence) and event=0 if censored: 
colon.surv <- na.omit(colon.dt)      
colon.surv[, event:=status.death]          
colon.surv[, time:=time.death]   
colon.surv <- colon.surv[, !(names(colon.surv) %in% grep("\\.death|\\.recurrence", names(colon.surv), value=TRUE)), with=FALSE]  

#-- create competing risks dataset with event=1 if cancer recurrence, event=2 if death without recurrence, event=0 if censored: 
colon.cr <- na.omit(colon.dt)                 
colon.cr[, time:=min(time.death, time.recurrence), by="id"]

fit.colon.fun <- function(formula.1=Surv(time, status==1)~rx+sex+nodes+age+obstruct+perfor+adhere+extent+surg,
                          formula.2=NULL,
                          formula.0=Surv(time, status==0)~rx+sex+nodes+age+obstruct+perfor+adhere+extent+surg,
                          formula.treat=trt~sex+age+nodes+obstruct+perfor+adhere+extent+surg,
                          initialize.scale.1=FALSE,
                          initialize.scale.2=FALSE,
                          initialize.scale.0=TRUE, #-- otherwise weibull for cens did not converge
                          d=colon) {

    covars <- unique(na.omit(gsub(" ", "", c(strsplit(as.character(formula.treat)[3], "\\+")[[1]],
                                             strsplit(as.character(formula.0)[3], "\\+")[[1]],
                                             strsplit(as.character(formula.1)[3], "\\+")[[1]],
                                             strsplit(as.character(formula.2)[3], "\\+")[[1]]))))

    covars.squared <- unlist(grep("\\.squared", covars))
    
    if (length(covars.squared)>0) {
        for (jj in covars.squared) {
            d[, (covars[jj]):=get(gsub("\\.squared", "", covars[jj]))^2]
        }
    }

    formula.treat.character <- as.character(formula.treat)
    d[, (paste0("num.", formula.treat.character[2])):=as.numeric(get(formula.treat.character[2]))]

    fit.treat <- fit.categorical(as.formula(paste0(paste0("num.", formula.treat.character[2]),
                                                   formula.treat.character[1],
                                                   formula.treat.character[3])), data=d, pseudo.out=TRUE)
    
    if (initialize.scale.1) {
        fit.weibull.1 <- survreg(as.formula(paste0(as.character(formula.1)[2], "~1")),
                                 data=na.omit(d), dist="weibull")
        fit.weibull.1 <- survreg(formula.1,
                                 data=na.omit(d), dist="weibull", scale=fit.weibull.1$scale)
    } else {
        fit.weibull.1 <- survreg(formula.1,
                                 data=na.omit(d), dist="weibull")
    }

    if (length(formula.2)>0) {
        if (initialize.scale.2) {
            fit.weibull.2 <- survreg(as.formula(paste0(as.character(formula.2)[2], "~1")),
                                     data=na.omit(d), dist="weibull")
            fit.weibull.2 <- survreg(formula.2,
                                     data=na.omit(d), dist="weibull", scale=fit.weibull.2$scale)
        } else {
            fit.weibull.2 <- survreg(formula.2,
                                     data=na.omit(d), dist="weibull")
        }
    }

    if (initialize.scale.0) {
        fit.weibull.0 <- survreg(as.formula(paste0(as.character(formula.0)[2], "~1")),
                                 data=na.omit(d), dist="weibull")
        fit.weibull.0 <- survreg(formula.0,
                                 data=na.omit(d), dist="weibull", scale=fit.weibull.0$scale)
    } else {
        fit.weibull.0 <- survreg(formula.0,
                                 data=na.omit(d), dist="weibull")
    }

    if (length(formula.2)>0) {
        return(list(fit.treat=fit.treat, fit.weibull.1=fit.weibull.1, fit.weibull.2=fit.weibull.2, fit.weibull.0=fit.weibull.0))
    } else {
        return(list(fit.treat=fit.treat, fit.weibull.1=fit.weibull.1, fit.weibull.0=fit.weibull.0))
    }

}

synthesize.colon.fun <- function(fit.colon, n=nrow(na.omit(d)), d=colon, get.true.value=NULL, tau=1000,
                                 name.treat="rx", event.name=NULL, time.name=NULL) {

    covars <- unique(unlist(lapply(fit.colon, function(fit) names(coef(fit))[-1])))
    covars <- covars[-unlist(lapply(c("pseudo", name.treat), function(x) grep(x, covars)))]

    d2 <- cbind(na.omit(d), fit.colon$fit.treat[[2]])
    
    covars2 <- names(d2)[unlist(sapply(names(d2), function(x) length(unique(grep(x, covars)>0))>0))]
    covars2 <- covars2[!covars2 %in% paste0(1:10)]
    
    sim.d <- setDT(d2)[sample(1:nrow(d2), n, replace=TRUE), c(covars2, colnames(fit.colon$fit.treat[[2]])), with=FALSE]
    
    if (length(get.true.value)>0) {
        sim.d[, (paste0(name.treat, ".num")):=get.true.value]
        sim.d[, (name.treat):=factor(get(paste0(name.treat, ".num")),
                                     levels=0:(length(levels(d2[[name.treat]]))-1),
                                     labels=levels(d2[[name.treat]]))]
    } else {
        sim.d[, (paste0(name.treat, ".num")):=as.numeric(sapply(1:nrow(sim.d), function(ii) sample(colnames(fit.colon$fit.treat[[2]]), 1, prob=sim.d[ii,names(sim.d)%in%colnames(fit.colon$fit.treat[[2]]), with=FALSE])))]
        sim.d[, (name.treat):=factor(get(paste0(name.treat, ".num")), labels=levels(d2[[name.treat]]))]
    }

    shape.1 <- 1/fit.colon$fit.weibull.1$scale
    scale.1 <- exp(predict(fit.colon$fit.weibull.1, newdata=sim.d, type="lp"))

    shape.0 <- 1/fit.colon$fit.weibull.0$scale
    scale.0 <- exp(predict(fit.colon$fit.weibull.0, newdata=sim.d, type="lp"))

    time.1 <- rweibull(n, shape=shape.1, scale=scale.1)
    time.0 <- rweibull(n, shape=shape.0, scale=scale.0)

    if ("fit.weibull.2" %in% names(fit.colon)) {
        shape.2 <- 1/fit.colon$fit.weibull.2$scale
        scale.2 <- exp(predict(fit.colon$fit.weibull.2, newdata=sim.d, type="lp"))
        time.2 <- rweibull(n, shape=shape.2, scale=scale.2)
        sim.d[, time:=sapply(1:n, function(i) min(time.2[i], time.1[i], time.0[i]))]
        sim.d[, status:=1*(time.1<=time.0 & time.1<=time.2) +
                    2*(time.2<=time.0 & time.2<time.1)]
        if (length(get.true.value)>0) {
            status.uncensored <- 1*(time.1<=time.2) + 2*(time.2<time.1)
            return(c(F1=mean(time.1<=tau & status.uncensored==1),
                     F2=mean(time.1<=tau & status.uncensored==2)))
        } 
    } else {
        sim.d[, time:=sapply(1:n, function(i) min(time.1[i], time.0[i]))]
        sim.d[, status:=1*(time.1<=time.0)]
        if (length(get.true.value)>0) {
            return(mean(time.1<=tau))
        } 
    }

    if (length(time.name)>0) {
        names(sim.d)[names(sim.d)=="time"] <- time.name
    }
    if (length(event.name)>0) {
        names(sim.d)[names(sim.d)=="status"] <- event.name
    }
        
    return(sim.d[, -colnames(fit.colon$fit.treat[[2]]), with=FALSE])
}

## Wrapper for simulating simpler stuff with binary treatment
colon_simulator <- function(formula_1,
                            formula_2,
                            formula_cens,
                            formula_treat){
    colon_binary <- copy(colon.cr)[, treat := factor(1*(rx != "Obs"))][, rx := NULL]
    fit.colon.cr <- fit.colon.fun(      
        formula.1=formula_1,
        formula.2=formula_2,
        formula.0=formula_cens,
        formula.treat=formula_treat,
        d=colon_binary)
    out = function(n, get_true_ate = FALSE, t = 1000){
        if(get_true_ate){
            true_A1 = as.numeric(synthesize.colon.fun(fit.colon=fit.colon.cr,
                                                      get.true.value = 1,
                                                      d=colon_binary,
                                                      name.treat="treat",
                                                      event.name="status",
                                                      n = n,
                                                      tau = t)["F1"])
            true_A0 = as.numeric(synthesize.colon.fun(fit.colon=fit.colon.cr,
                                                      get.true.value = 0,
                                                      d=colon_binary,
                                                      name.treat="treat",
                                                      event.name="status",
                                                      n = n,
                                                      tau = t)["F1"])
            true_ate <- true_A1-true_A0
            return(list("A=1" = true_A1, "A=0" = true_A0, "ATE" = true_ate))
        }else{
            sim_dat0 <- synthesize.colon.fun(fit.colon=fit.colon.cr,
                                             d=colon_binary,
                                             name.treat="treat",
                                             event.name="status",
                                             n = n)
            ## Setup data on right format
            sim_dat0[, A := as.numeric(treat)-1][, ":="(treat.num = NULL, treat = NULL)]
            return(sim_dat0[])
        }
    }
    return(out)
}


######################################################################
### hely-colon-sim.R ends here
