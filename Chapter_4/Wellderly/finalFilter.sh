#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

module load vcftools

vcftools --vcf wellderlyAll.recode.vcf --keep indivToKeep.txt --recode --out wellderlyFiltered

