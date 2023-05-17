#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 750
#SBATCH --constraint=rhel8
#SBATCH --output=./../cluster-out/equal-samp-log-gamma-re-%a.out
#SBATCH --error=./../cluster-err/equal-samp-log-gamma-re-%a.err
#SBATCH --array=1-2500

## add R module
module add r/4.1.3

## install R package - execute the following line only once
## Rscript --vanilla ./../R-scripts/install_packages.R

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-equal-samp-log-gamma-re.R ./../cluster-log/equal-samp-log-gamma-re/equal-samp-log-gamma-re-$SLURM_ARRAY_TASK_ID.Rout
