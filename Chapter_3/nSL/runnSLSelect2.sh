#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

for i in {1..100}
do
	nSL -adfile nSLSelectedInputFiles2/adfile_$i.txt -samfile nSLSelectedInputFiles2/samfile_$i.txt -hapfile nSLSelectedInputFiles2/hapfile_$i.txt > nSLSelectedOutputFiles2/nslSelect_$i.txt
done
