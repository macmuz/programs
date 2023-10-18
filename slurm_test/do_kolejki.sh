#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --ntasks-per-node=24
for i in `seq 24`; do
  srun --exclusive --ntasks 1 ./test ${i} &
done
wait
