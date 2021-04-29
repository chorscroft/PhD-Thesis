#!/bin/bash
#PBS -l walltime=00:02:00
cd $PBS_O_WORKDIR

module load R/3.2.1
Rscript createAggregateGraphsNew.R
