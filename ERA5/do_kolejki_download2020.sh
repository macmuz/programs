#!/bin/bash -f

#SBATCH -J ERA
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=16:00:00

python get_year_ERA5_baltic.py 2020