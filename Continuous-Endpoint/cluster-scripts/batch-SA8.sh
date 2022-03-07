#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 1000
#SBATCH --output=./../cluster-out/SA8-%a.out
#SBATCH --error=./../cluster-err/SA8-%a.err
#SBATCH --array=1-60

## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-SA8.R ./../cluster-log/SA8/SA8-$SLURM_ARRAY_TASK_ID.Rout