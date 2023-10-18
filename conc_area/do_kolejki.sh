#!/bin/bash -f

#SBATCH -J conc
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00

./conc_area
