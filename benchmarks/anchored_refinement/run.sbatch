#!/bin/bash
#SBATCH --job-name=mr4-anchored
#SBATCH --array=1-200
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/rep_%a.log

## Slurm launcher for REG-MR4-anchored-v1 (staged; NOT launched from the
## sandbox). Requires: cellprograms + flashier.utils + clue + ebnm on the
## module path; the repo checked out at $CP_REPO; amplitude pilot completed
## and MR4_AMP exported (or amplitudes.csv sourced below).

export CP_REPO=${CP_REPO:-/scratch/cellprograms}
mkdir -p logs out
export MR4_AMP=${MR4_AMP:-1}

Rscript run_rep.R $SLURM_ARRAY_TASK_ID out
