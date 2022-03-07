#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 3000
#SBATCH --output=./../cluster-out/n7650-beta20-%a.out
#SBATCH --error=./../cluster-err/n7650-beta20-%a.err
#SBATCH --array=1-1200

## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-n7650-beta20.R ./../cluster-log/n7650-beta20/n7650-beta20-$SLURM_ARRAY_TASK_ID.Rout