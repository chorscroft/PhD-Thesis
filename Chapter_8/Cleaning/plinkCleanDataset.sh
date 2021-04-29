#!/bin/bash
#PBS -l walltime=00:02:00
cd $PBS_O_WORKDIR


## Cleans the All_3_merged dog dataset

module load plink/1.90beta
plink --dog --file /temp/hgig/EXOME_DATA/LD_MAPS/Dog/merge/raw_data/All_3_merged_raw --maf 0.05 --geno 0.05 --mind 0.05 --snps-only no-DI --extract /temp/hgig/EXOME_DATA/LD_MAPS/Dog/SOURCE_DATA/sams_boyko_ROH_2018_data/IDs --make-bed --out /temp/hgig/EXOME_DATA/Clare/Dogs/CleanDataset/clean
plink --dog --hardy --bfile /temp/hgig/EXOME_DATA/Clare/Dogs/CleanDataset/clean
