#!/bin/bash
#PBS -l walltime=07:00:00

cd $PBS_O_WORKDIR

module load jdk
java -jar msms.jar 100 100 -t 2000 -r 19999.96 500000 > msmsNeutral.txt

