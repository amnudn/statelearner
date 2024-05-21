library(data.table)
options(rf.cores = 1) ## Keep RF on one core
setDTthreads(1) ## Avoid DT spawning

## Get input
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
max_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_MAX"))
job_name <- as.character(Sys.getenv("SLURM_JOB_NAME"))

## Get seeds
set.seed(284219)
seeds <- sample(1:1e6, size = max_task_id, replace = FALSE)

## Load libraries here:
library(here)
library(targets)
tar_source(here("R-code/functions"))

## Sample one draw from each estimator and simulation setting
result <- numerical_study_sampler(seeds[task_id])

## Save results externally to results folder
fwrite(result, file = paste0('./results/', job_name, '_taskID_', sprintf("%05d", task_id), '.txt'))
