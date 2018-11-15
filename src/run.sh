#!/usr/bin/env bash

#SBATCH -J paper
#SBATCH -p DPB
#SBATCH -c 2
#SBATCH --mem=4000

module load Java/8 Trinity/2.3.2 Bowtie/2.2.9 Python/3.6.0

srun python3 paper.py -a -i $1 -o $2 -c 12 -m 95
