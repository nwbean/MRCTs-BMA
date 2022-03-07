#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 1000
#SBATCH --output=./../cluster-out/SA9-%a.out
#SBATCH --error=./../cluster-err/SA9-%a.err
#SBATCH --array=1-60

## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-SA9.R ./../cluster-log/SA9/SA9-$SLURM_ARRAY_TASK_ID.Rout