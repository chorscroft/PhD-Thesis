

#!/bin/bash

#PBS -l walltime=00:10:00

cd $PBS_O_WORKDIR


/home/ch19g17/LDHat/LDHatMaster/stat -input rates.txt -loc locs.txt -burn 20
