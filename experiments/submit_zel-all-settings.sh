#!/bin/bash

#SBATCH --job-name=zel-all-settings
#SBATCH --array=1-300%60
#SBATCH --output=log/stdout-%x-%a.out

# Enable Additional Software
module load gcc/11.2.0
module load R/4.1.2

# Run the job commands
Rscript --vanilla $SLURM_JOB_NAME.R
