#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

for i in {1..100}
do
	nSL -adfile nSLSelectedInputFiles/adfile_$i.txt -samfile nSLSelectedInputFiles/samfile_$i.txt -hapfile nSLSelectedInputFiles/hapfile_$i.txt > nSLSelectedOutputFiles/nslSelect_$i.txt
done
