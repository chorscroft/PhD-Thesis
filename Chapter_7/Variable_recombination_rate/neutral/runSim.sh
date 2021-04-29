#!/bin/bash
#PBS -l walltime=00:10:00

cd $PBS_O_WORKDIR

/home/ch19g17/Sim_for_variable_recombination_zalpha/build/slim neutralsim.txt
