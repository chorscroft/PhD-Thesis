#!/bin/bash
#PBS -l walltime=00:30:00
cd $PBS_O_WORKDIR

module load vcftools
vcftools --vcf /temp/hgig/EXOME_DATA/Clare/wellderlyCh2LCT.recode.vcf --min-alleles 2 --max-alleles 2 --plink-tped --out wellderlyChr2LCTtped

#vcftools --vcf /temp/hgig/EXOME_DATA/Clare/wellderlyCh2LCT.recode.vcf --chr chr2 --from-bp 130000000 --to-bp 130050000 --min-alleles 2 --max-alleles 2 --plink-tped --out wellderlyChr2LCTtpedSmall
