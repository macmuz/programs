#!/bin/bash -f

#SBATCH -J ROMS_bc
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

./ROMS_bc_spline 2016 ROMS_grid_05NM_015p3.nc
