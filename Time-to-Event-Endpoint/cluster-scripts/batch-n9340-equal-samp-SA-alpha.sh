#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 750
#SBATCH --output=./../cluster-out/n9340-equal-samp-SA-alpha-%a.out
#SBATCH --error=./../cluster-err/n9340-equal-samp-SA-alpha-%a.err
#SBATCH --array=1-100

## add R module
module add r/4.0.3

## install R package - execute the following line only once
## Rscript --vanilla ./../R-scripts/install_packages.R

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-n9340-equal-samp-SA-alpha.R ./../cluster-log/n9340-equal-samp-SA-alpha/n9340-equal-samp-SA-alpha-$SLURM_ARRAY_TASK_ID.Rout
