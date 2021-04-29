
#!/bin/bash

#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR


module load plink/1.90beta

plink --file wellderlyAllpruned --pca var-wts
