

#!/bin/bash

#PBS -l walltime=00:10:00

cd $PBS_O_WORKDIR


../LDHatMaster/stat -input rates.txt -loc CEU22.locs -burn 20
