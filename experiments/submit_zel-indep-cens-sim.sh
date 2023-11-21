#!/bin/bash

#SBATCH --job-name=zel-indep-cens-sim
#SBATCH --array=1-500%40
#SBATCH --output=log/stdout-%x-%a.out

# Enable Additional Software
module load gcc/11.2.0
module load R/4.1.2

# Run the job commands
Rscript --vanilla zel-indep-cens-sim.R
