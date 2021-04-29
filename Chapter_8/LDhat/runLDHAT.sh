#!/bin/bash

#PBS -l walltime=02:00:00

cd $PBS_O_WORKDIR

chromDir=/temp/hgig/EXOME_DATA/Clare/Dogs/LDprofile/LDHAT/chr$chrom

## Converts the clipped, clean dog data into VCF
module load plink/1.90beta
plink --dog --tfile /temp/hgig/EXOME_DATA/Clare/Dogs/zalpha/DataForZalpha/chr$chrom --chr $chrom --recode vcf --out $chromDir/cleanVCF

## Creates the input files for LDhat
module load R/3.6.1

Rscript $chromDir/convertVcfDataToLDHatInputFiles.R $chromDir --vanilla

## Run the convert module of LDhat
/home/ch19g17/LDHat/LDHatMaster/convert -seq $chromDir/dogs.sites -loc $chromDir/dogs.locs -nout 96

## Run the interval module of LDhat
/home/ch19g17/LDHat/LDHatMaster/interval -seq $chromDir/sites.txt -loc $chromDir/locs.txt -lk /home/ch19g17/LDHat/lk_n192_t0.001 -its 10000000 -bpen 5 -samp 5000

## Run the stat module of LDhat
/home/ch19g17/LDHat/LDHatMaster/stat -input $chromDir/rates.txt -loc $chromDir/locs.txt -burn 20
