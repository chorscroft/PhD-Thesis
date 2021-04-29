#!/bin/bash

#PBS -l walltime=01:00:00
cd $PBS_O_WORKDIR

module load R/3.6.1

R CMD BATCH findRelatedDogsToRemove.R
