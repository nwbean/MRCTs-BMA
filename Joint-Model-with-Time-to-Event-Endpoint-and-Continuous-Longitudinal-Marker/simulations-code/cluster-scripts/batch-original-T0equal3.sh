#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 750
#SBATCH --constraint=rhel8
#SBATCH --output=./../cluster-out/original-T0equal3-%a.out
#SBATCH --error=./../cluster-err/original-T0equal3-%a.err
#SBATCH --array=1-4000

## add R module
module add r/4.1.3

## install R package - execute the following line only once
## Rscript --vanilla ./../R-scripts/install_packages.R

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" ./../R-scripts/run-simulation-original-T0equal3.R ./../cluster-log/original-T0equal3/original-T0equal3-$SLURM_ARRAY_TASK_ID.Rout
