library(riskRegression)
library(data.table)
library(here)

set.seed(12313)
emulated_zel <- simZelefsky(n = 1042)[, .(time = dmos,status = 1*status,logPSA,stage,ggtot,sDose,hormones)]
emulated_zel[, death := rexp(n = .N, rate = 1/500)]
emulated_zel[death <= time, status := 2]
emulated_zel[death <= time, time := death]
emulated_zel[, death := NULL]

fwrite(emulated_zel, file = here("emulated-data/emulated-data.txt"))
