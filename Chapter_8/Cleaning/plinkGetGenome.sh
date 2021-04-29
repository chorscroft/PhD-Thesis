#!/bin/bash

#PBS -l walltime=00:20:00
cd $PBS_O_WORKDIR

module load plink/1.90beta

plink --bfile ../Prune/pruned --dog --genome
