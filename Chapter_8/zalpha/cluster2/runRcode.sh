#!/bin/bash

#PBS -l walltime=02:00:00
cd $PBS_O_WORKDIR

module load R/3.6.1

Rscript runZalpha.R $j --vanilla
