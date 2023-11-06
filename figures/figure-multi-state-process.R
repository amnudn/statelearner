### figure-multi-state-process.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 27 2023 (10:04) 
## Version: 
## Last-Updated: Oct  9 2023 (18:22) 
##           By: Anders Munch
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(prodlim)
nTrans <- 3
stateLabels = c("Initial","Cause 1", "Cause 2", "Censored")
crHist <- Hist(time = 1:nTrans, event = list(from = rep("1", nTrans), to = stateLabels[-1]))
pdf("~/research/SuperVision/Anders/survival-loss/statelearner/figure-multi-state-process.pdf")
plot(crHist,stateLabels = c("Initial","Cause 1", "Cause 2", "Censored"),arrowLabels = FALSE)
dev.off()

######################################################################
### figure-multi-state-process.R ends here
