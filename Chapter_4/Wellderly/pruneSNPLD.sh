
#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR


module load plink/1.90beta

plink --vcf wellderlyAll.recode.vcf --indep-pairwise 50 5 0.5

plink --vcf wellderlyAll.recode.vcf --extract plink.prune.in --recode --out wellderlyAllpruned
