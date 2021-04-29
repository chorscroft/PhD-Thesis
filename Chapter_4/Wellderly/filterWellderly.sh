#!/bin/bash
#PBS -l walltime=05:00:00
cd $PBS_O_WORKDIR
module load vcftools
vcftools --vcf /temp/hgig/EXOME_DATA/LD_MAPS/CG_Wellderly/rawdata/Wellderly_597-Small_Variant_Table-SNPs.tsv.vcf --max-missing 0.9 --maf 0.05 --recode --hwe 0.001 --out wellderlyAll


