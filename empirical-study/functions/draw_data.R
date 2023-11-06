### draw_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2023 (13:07) 
## Version: 
## Last-Updated: Jul 11 2023 (14:15) 
##           By: Thomas Alexander Gerds
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
draw_data <- function(n){
    ## set effects
    beta1 <- rep(c(1,-1,0,0,1,-1),2) * log(2)
    q.beta1 <- c(1,-1,0,0,1,-1) * log(2)
    q.gamma <- c(0,0,1,-1,1,-1) * log(2)
    gamma <- rep(c(0,0,1,-1,1,-1),2) * log(2)
    ## create data generating model 
    qmodel <- quadratic_model(beta1=beta1,q.vars=paste("X",1:6,sep=""),q.beta1=q.beta1,q.gamma=q.gamma,gamma=gamma,competing.risks = FALSE)
    d = sim(qmodel,n)
    # add censoring status
    setDT(d)
    d[,cens_status := 1*(event != 0)]
    d[]
}


######################################################################
### draw_data.R ends here
