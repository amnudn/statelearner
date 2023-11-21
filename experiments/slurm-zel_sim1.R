library(data.table)
options(rf.cores = 1) ## Keep RF on one core
setDTthreads(1) ## Avoid DT spawning

## Get input
number_of_tasks <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Get seeds
set.seed(12324)
seeds <- sample(1:1e6, size = 1e4, replace = FALSE)

## Run sim function with current seed
source("zel-based-sim.R")
result <- zel_sim1(seeds[task_id])

# Save the results for this task as an individual file in the output folder
## if(!dir.exists("./results")) dir.create("./results")
fwrite(result, file = paste0('./results/zel_sim1-', sprintf("%05d", task_id), '.txt'))

## TODO: make it end with -task-id
