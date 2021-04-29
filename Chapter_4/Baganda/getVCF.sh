
#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR


module load plink/1.90beta

plink --bfile chr22_Baganda --maf 0.05 --hwe 0.001 --geno 0.1 --biallelic-only --recode vcf --out chr22_Baganda
