#!/bin/bash

#PBS -l walltime=60:00:00

cd $PBS_O_WORKDIR

../LDHatMaster/interval -seq CEU22.sites -loc CEU22.locs -lk ../lk_n192_t0.001 -its 10000000 -bpen 5 -samp 5000
