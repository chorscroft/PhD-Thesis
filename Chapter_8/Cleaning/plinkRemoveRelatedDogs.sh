#!/bin/bash

#PBS -l walltime=00:02:00
cd $PBS_O_WORKDIR

module load plink/1.90beta

plink --bfile ../Prune/pruned --dog --remove dogstoremove.txt --make-bed --out unrelated
