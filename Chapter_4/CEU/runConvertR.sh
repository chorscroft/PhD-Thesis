#!/bin/bash

#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

module load R/3.4.2

R --file=convertVcfDataToLDHatInputFiles.R
#R --file=convertMultipleFiles.R

