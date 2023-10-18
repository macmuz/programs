#!/bin/bash -f

#SBATCH -J diag
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

./section ocean_avg_2004-02-16.nc
