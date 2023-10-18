#!/bin/bash -f

#SBATCH -J ROMS_force
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00

#./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_Tair_2019.nc tair_time Tair &
./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_wind_2019.nc wind_time Uwind Vwind
#./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_lwrad_down_2019.nc lrf_time lwrad_down &
#./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_Pair_2019.nc pair_time Pair &
#./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_Qair_2019.nc qair_time Qair &
#./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_rain_2019.nc rain_time rain &
#./move_forcing /users/work/mmuzyka/CSDIR/forcing/baltic_swrad_2019.nc srf_time swrad &

#wait
