### figure-multi-state-process.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 27 2023 (10:04) 
## Version: 
## Last-Updated: Apr 29 2024 (10:30) 
##           By: Anders Munch
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(prodlim)
library(here)
nTrans <- 3
stateLabels = c("Initial","Cause 1", "Cause 2", "Censored")
crHist <- Hist(time = 1:nTrans, event = list(from = rep("1", nTrans), to = stateLabels[-1]))
try(setwd("~/research/SuperVision/Anders/survival-loss/statelearner/figures"))
try(setwd(here("figures")))
pdf("figure-multi-state-process.pdf")
plot(crHist,stateLabels = stateLabels,arrowLabels = FALSE,
     tagBoxes = FALSE) ## Make true and/or change to -1?
dev.off()

######################################################################
### figure-multi-state-process.R ends here
