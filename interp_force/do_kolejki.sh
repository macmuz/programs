#!/bin/bash

#SBATCH -J interp
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

cd /users/work/mmuzyka/programs/interp_force

./calc lwrad.in &
./calc pair.in &
./calc qair.in &
./calc rain.in &
./calc swrad.in &
./calc tair.in &
./calc wind.in &
wait
