#!/bin/bash
#PBS -l walltime=05:00:00
cd $PBS_O_WORKDIR
module load vcftools
vcftools --vcf /temp/hgig/EXOME_DATA/LD_MAPS/CG_Wellderly/rawdata/Wellderly_597-Small_Variant_Table-SNPs.tsv.vcf --maf 0.01 --chr chr2 --recode --hwe 0.001 --out wellderlyCh2


