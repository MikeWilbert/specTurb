#!/bin/bash -x
#SBATCH --account=specturb
#SBATCH --nodes=6
#SBATCH --tasks-per-node=48
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --time=00:10:00
#SBATCH --partition=devel

srun --overlap ./TSpecDyn
