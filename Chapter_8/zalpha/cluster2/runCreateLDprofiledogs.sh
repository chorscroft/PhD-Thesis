#!/bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=12

cd $PBS_O_WORKDIR

module load R/3.6.1
Rscript Create_LDprofile_Dogs.R --vanilla
