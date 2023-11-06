### minus_t_fun.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  6 2023 (12:11) 
## Version: 
## Last-Updated: Nov  6 2023 (12:13) 
##           By: Anders Munch
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary:
##
## From
##    t \mapsto Lambda(newdata, t),
## construct the function
##    t \mapsto Lambda(newdata, t-).
##
## It is assumed that Lambda only jumps at jump_points.
## 
## Returns a function: out=function(newdata,time)
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## Meta functions:
minus_t_fun <- function(Lambda, jump_points){
    t_minus = min(diff(sort(unique(jump_points))))/2
    fun0 = function(newdata, times){
        times0 = pmax(0, times-t_minus)
        Lambda(newdata, times0)
    }
    return(fun0)
}


######################################################################
### minus_t_fun.R ends here
