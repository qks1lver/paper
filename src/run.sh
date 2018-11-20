#!/usr/bin/env bash

#SBATCH -J paper
#SBATCH -p DPB
#SBATCH -c 2
#SBATCH --mem=4000

module load Python/3.6.0

srun python3 paper.py "$@" 2>&1
