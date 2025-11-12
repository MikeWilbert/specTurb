#!/bin/bash -x
#SBATCH --account=specturb
#SBATCH --nodes=192
#SBATCH --tasks-per-node=48
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=batch

srun --overlap ./TSpecDyn_B10
