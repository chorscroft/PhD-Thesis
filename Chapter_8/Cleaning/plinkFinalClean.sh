#!/bin/bash
#PBS -l walltime=00:02:00
cd $PBS_O_WORKDIR


## Creates the final dataset

module load plink/1.90beta
plink --dog --bfile /temp/hgig/EXOME_DATA/Clare/Dogs/CleanDataset/clean --keep /temp/hgig/EXOME_DATA/Clare/Dogs/CleanDataset/Relatedness/PCA/indivToKeep.txt --maf 0.05 --geno 0.05 --recode --transpose --out cleanDogDataset
