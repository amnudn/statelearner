### quadratic_model.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2023 (13:08) 
## Version: 
## Last-Updated: Jul 11 2023 (13:18) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
quadratic_model <- function(beta1,
                            beta2,
                            gamma,
                            q.vars,
                            q.beta1,
                            q.beta2,
                            q.gamma,
                            competing.risks = TRUE){
    require(lava)
    m <- lvm()
    X <- paste("X",seq(12),sep="")
    timecause1 <- "timecause1"
    if (competing.risks[[1]]) timecause2 <- "timecause2"
    censtime <- "censtime"
    if (competing.risks[[1]])
        latent(m) <- c(timecause1,timecause2,censtime)
    else
        latent(m) <- c(timecause1,censtime)
    if (competing.risks[[1]])
        addvar(m) <- c(timecause1,timecause2,censtime,X)
    else
        addvar(m) <- c(timecause1,censtime,X)
    exogenous(m) <- c(X,q.vars)
    distribution(m,X[1:6]) <- normal.lvm()
    distribution(m,X[7:12]) <- binomial.lvm()
    distribution(m,timecause1) <- coxWeibull.lvm(scale=1/100)
    if (competing.risks[[1]]) distribution(m,timecause2) <- coxWeibull.lvm(scale=1/100)
    distribution(m,censtime) <- coxWeibull.lvm(scale=1/100)
    regression(m,to=timecause1,from=X) <- beta1
    if (competing.risks[[1]]) regression(m,to=timecause2,from=X) <- beta2
    regression(m,to=censtime,from=X) <- gamma
    ## quadratic terms
    for (v in q.vars){
        form.v <- formula(paste("quad.",v,"~",v,sep=""))
        transform(m,form.v) <- function(x)x^2
    }
    regression(m,to=timecause1,from=paste("quad.",q.vars,sep="")) <- q.beta1
    if (competing.risks[[1]]) regression(m,to=timecause2,from=paste("quad.",q.vars,sep="")) <- q.beta2
    regression(m,to=censtime,from=paste("quad.",q.vars,sep="")) <- q.gamma
    if (competing.risks[[1]])
        m <- eventTime(m,time~min(timecause1=1,timecause2=2,censtime=0),"event")
    else
        m <- eventTime(m,time~min(timecause1=1,censtime=0),"event")
    return(m)
}
######################################################################
### quadratic_model.R ends here
