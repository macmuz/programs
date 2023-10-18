#!/bin/bash -f

#SBATCH -J ROMS_bc
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=16:00:00

./new_ROMS_bc_bar 2023 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc 
