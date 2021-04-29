#!/bin/bash

#PBS -l walltime=13:00:00

cd $PBS_O_WORKDIR

module load R/3.4.2

R --file=analyseH12SimP3.R -q
