#!/bin/bash

#PBS -l walltime=01:00:00

cd $PBS_O_WORKDIR

## start loop
for i in {1..38}  #38
do
	## Make directory	
	mkdir chr$i

	## Copy code
	cp runLDHAT.sh chr$i/
	cp convertVcfDataToLDHatInputFiles.R chr$i/
	
	cd chr$i/
	## set code running
	qsub runLDHAT.sh -v chrom=$i
	cd ..

done



