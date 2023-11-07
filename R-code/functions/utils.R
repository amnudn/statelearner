### utils.R --- 
#----------------------------------------------------------------------
## Author: Anders Munch
## Created: Nov  7 2023 (10:02) 
## Version: 
## Last-Updated: Nov  7 2023 (10:02) 
##           By: Anders Munch
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## util for rewrapping list of lists 
rewrap <- function(ll, comb_fun = rbind){
    elements_in_list = length(ll[[1]])
    names_of_elements = names(ll[[1]])
    if(class(comb_fun) == "function")
        rep(c(comb_fun), times = elements_in_list)
    stopifnot(elements_in_list == length(comb_fun))
    out = lapply(1:elements_in_list, function(ii){
        do.call(comb_fun[[ii]], lapply(ll, function(xx){
            xx[[ii]]
        }))
    })
    names(out) = names_of_elements
    return(out)
}



######################################################################
### utils.R ends here
