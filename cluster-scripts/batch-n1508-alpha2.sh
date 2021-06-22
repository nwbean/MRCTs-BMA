#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 1000
#SBATCH --output=./../cluster-out/n1508-alpha2-%a.out
#SBATCH --error=./../cluster-err/n1508-alpha2-%a.err
#SBATCH --array=1-60

## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-n1508-alpha2.R ./../cluster-log/n1508-alpha2/n1508-alpha2-$SLURM_ARRAY_TASK_ID.Rout