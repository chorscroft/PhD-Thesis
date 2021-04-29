## Runs the LDHat code to "convert"
## chooses 96 sequences at random

#!/bin/bash

#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR

/home/ch19g17/LDHat/LDHatMaster/convert -seq Baganda22.sites -loc Baganda22.locs -2only -freqcut 0.05 -missfreqcut 0.1 -nout 96
