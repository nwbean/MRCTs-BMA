#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 1000
#SBATCH --output=./../cluster-out/SA7-%a.out
#SBATCH --error=./../cluster-err/SA7-%a.err
#SBATCH --array=1-60

## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-SA7.R ./../cluster-log/SA7/SA7-$SLURM_ARRAY_TASK_ID.Rout