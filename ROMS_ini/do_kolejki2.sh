#!/bin/bash -f

#SBATCH -J ROMS_ini
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

./ROMS_ini ROMS_grid_125NM_new_hmask_25.nc 1992 7 1

