### at_risk_fun.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  6 2023 (12:07) 
## Version: 
## Last-Updated: Nov  6 2023 (12:10) 
##           By: Anders Munch
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary:
##
## In newdata, use time_name to calculate who are at risk at
## each time point in time.
##
## Returns matrix of dimension nrow(newdata) * length(times)
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)

at_risk_fun <- function(newdata,times,time_name = "time"){
    spring_times = newdata[[time_name]]
    grid_mat = matrix(times,ncol = length(times),nrow = length(spring_times),byrow = TRUE)
    ## event_mat = apply(grid_mat, 2, function(x) 1*(spring_times >= x))
    event_mat = do.call(cbind, lapply(1:ncol(grid_mat), function(xx) 1*(spring_times >= grid_mat[, xx])))
    return(event_mat)
}



######################################################################
### at_risk_fun.R ends here
