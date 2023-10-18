#!/bin/bash -f

#SBATCH -J ROMS_ini
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

#./ROMS_ini /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc 2003 11 1
./ROMS_ini /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_2_3km_560x600_NetCDF4.nc 1993 1 2

