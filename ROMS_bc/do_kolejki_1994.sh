#!/bin/bash -f

#SBATCH -J ROMS_bc
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

./old_ROMS_bc_bar 1994 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_2_3km_560x600_NetCDF4.nc 
