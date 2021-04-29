#!/bin/bash
#PBS -l walltime=14:00:00

cd $PBS_O_WORKDIR

module load jdk
java -jar msms.jar 100 100 -t 2000 -r 19999.96 500000 -N 1000000 -SAA 24000 -SAa 12000 -SF 0 0.5 -Sp 0.5 -seed 0x6c0db7a08ee8a468 -T  > msmsSelectTreePT.txt

