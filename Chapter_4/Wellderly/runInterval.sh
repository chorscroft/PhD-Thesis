#!/bin/bash

#PBS -l walltime=20:00:00

cd $PBS_O_WORKDIR

/home/ch19g17/LDHat/LDHatMaster/interval -seq sites.txt -loc locs.txt -lk /home/ch19g17/LDHat/lk_n192_t0.001 -its 10000000 -bpen 5 -samp 5000
