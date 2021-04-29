#!/bin/bash

#PBS -l walltime=00:10:00
cd $PBS_O_WORKDIR

module load plink/1.90beta

## get list of boundaries to filter by from Auton map limits
## output .tped files for each chromosome

i=1
while read -ra f
do
	## filter by them and create new .tpeds for each chromosome
	plink --tfile /temp/hgig/EXOME_DATA/Clare/Dogs/CleanDataset/FinalCleanDataset/cleanDogDataset --dog --chr $i --from-bp ${f[1]} --to-bp ${f[2]} --recode transpose --out chr$i
	
	##iterate chromosome
	i=$((i+1))

done < /temp/hgig/EXOME_DATA/Clare/Dogs/zalpha/DataForZalpha/AutonBoundaries.txt


