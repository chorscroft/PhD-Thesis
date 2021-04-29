
#!/bin/bash

#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR
module load vcftools

vcftools --vcf ../WellderlyRemoveOutliers/wellderlyFiltered.recode.vcf --chr chr22 --max-alleles 2 --recode --out wellderlyChr22
 
