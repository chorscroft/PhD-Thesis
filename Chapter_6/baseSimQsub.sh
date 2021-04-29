#!/bin/bash
#PBS -l walltime=02:00:00

cd $PBS_O_WORKDIR

../build/slim SIMTORUN.txt

