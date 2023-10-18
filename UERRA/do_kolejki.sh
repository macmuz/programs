#!/bin/bash -f

#SBATCH -J ROMS_bc
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00

YEAR=2011

./cds_request_swrad_down.py $YEAR &
./cds_request_lwrad_down.py $YEAR &
./cds_request_swrad.py $YEAR &
./cds_request_rain.py $YEAR &
./cds_request_wdir.py $YEAR &
./cds_request_wspeed.py $YEAR &
./cds_request_tair.py $YEAR &
./cds_request_pair.py $YEAR &
./cds_request_qair.py $YEAR &
wait
