#!/bin/bash -f

#SBATCH -J ROMS_bc
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:15:00

./old_ROMS_bc_bar 2014 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2015 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2016 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2017 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2018 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2019 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2020 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
./old_ROMS_bc_bar 2021 /users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc &
wait
