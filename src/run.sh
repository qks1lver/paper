#!/usr/bin/env bash

#SBATCH -J paper
#SBATCH -p DPB
#SBATCH -c 1
#SBATCH --mem=2000

srun python3 paper.py -a -i $1 -o $2 -c $3 -m $4
