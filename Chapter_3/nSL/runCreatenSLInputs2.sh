
#!/bin/bash

#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR


module load R/3.4.2

R --file=startCreatingInputs2.R -q


