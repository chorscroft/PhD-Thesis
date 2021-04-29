#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

for i in {1..100}
do
	nSL -adfile nSLSelectedInputFilesP3/adfile_$i.txt -samfile nSLSelectedInputFilesP3/samfile_$i.txt -hapfile nSLSelectedInputFilesP3/hapfile_$i.txt > nSLSelectedOutputFilesP3/nslSelect_$i.txt
done
