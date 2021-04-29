#!/bin/bash
#PBS -l walltime=09:00:00
cd $PBS_O_WORKDIR

module load R/3.4.2
Rscript getResults_wideROC.R
