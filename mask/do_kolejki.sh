#!/bin/bash -f

#SBATCH -J mask
#SBATCH --partition batch_16h
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=16:00:00

./mask
