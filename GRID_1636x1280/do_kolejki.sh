#!/bin/bash -f

#SBATCH -J HRDM_bathy
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00

./HRDM_bathy >log
