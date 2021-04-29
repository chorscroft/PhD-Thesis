
#!/bin/bash

#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR


module load plink/1.90beta

plink --bfile ../IBDgenome/unrelated --pca 3701 --make-rel --dog --chr 1-38

awk '{print $NR}' plink.rel > plink.rel.diag
