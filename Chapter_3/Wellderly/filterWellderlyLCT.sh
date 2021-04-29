
#!/bin/bash
#PBS -l walltime=02:00:00
cd $PBS_O_WORKDIR

#This code filters the file by bp 135000000 to 138000000 on Chr2

module load vcftools

vcftools --vcf wellderlyCh2.recode.vcf --chr chr2 --from-bp 135000000 --to-bp 138000000 --recode --out wellderlyCh2LCT
