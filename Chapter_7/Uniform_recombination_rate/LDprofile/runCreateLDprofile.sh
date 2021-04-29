#!/bin/bash
#PBS -l walltime=00:15:00
cd $PBS_O_WORKDIR

module load R/3.4.2
Rscript createLDprofile.R
