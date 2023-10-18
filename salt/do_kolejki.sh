#!/bin/bash -f

#SBATCH -J salt
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=10:00:00

./salt_ROMS
