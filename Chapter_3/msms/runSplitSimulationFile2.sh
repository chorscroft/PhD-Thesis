#!/bin/bash

#PBS -l walltime=00:07:00

cd $PBS_O_WORKDIR

module load R/3.4.2

R --file=splitSimulationFile2.R -q

