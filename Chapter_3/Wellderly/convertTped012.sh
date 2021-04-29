#!/bin/bash

#PBS -l walltime=00:01:00
cd $PBS_O_WORKDIR

module load plink/1.90beta

plink --tfile wellderlyChr2LCTtped --recode12 --transpose --out wellderly012
