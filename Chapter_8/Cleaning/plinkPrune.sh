#!/bin/bash

#PBS -l walltime=00:20:00
cd $PBS_O_WORKDIR

module load plink/1.90beta

plink --bfile ../../clean --indep-pairwise 50 5 0.5 --dog --out pruned
plink --bfile ../../clean --extract pruned.prune.in --dog --make-bed --out pruned
