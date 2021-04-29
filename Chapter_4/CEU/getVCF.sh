
#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

module load vcftools

vcftools --vcf Chr22_CEU_nohemiRef_20mb+_maf0.05-0.95_miss0.95_hwe0.01.recode.vcf --max-alleles 2 --recode --out CEU22

 
