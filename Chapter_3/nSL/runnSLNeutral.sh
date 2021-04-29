#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

for i in {1..100}
do
	nSL -adfile nSLNeutralInputFiles/adfile_$i.txt -samfile nSLNeutralInputFiles/samfile_$i.txt -hapfile nSLNeutralInputFiles/hapfile_$i.txt > nSLNeutralOutputFiles/nslNeutral_$i.txt
done
