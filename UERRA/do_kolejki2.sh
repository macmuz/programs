#!/bin/bash -f

#SBATCH -J UERRA
#SBATCH -p batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00

./forcing 2008
