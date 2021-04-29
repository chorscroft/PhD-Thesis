#!/bin/bash

#PBS -l walltime=05:00:00

cd $PBS_O_WORKDIR

for i in {1..100}
do
	nSL -adfile nSLSelectedInputFilesP/adfile_$i.txt -samfile nSLSelectedInputFilesP/samfile_$i.txt -hapfile nSLSelectedInputFilesP/hapfile_$i.txt > nSLSelectedOutputFilesP/nslSelect_$i.txt
done
