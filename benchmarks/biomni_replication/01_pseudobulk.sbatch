#!/bin/bash
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --job-name=pb
#SBATCH --output=results/logs/pb_%x_%j.log
# usage: sbatch [-J pb_<cohort>] 01_pseudobulk.sbatch <cohort> [--source X|raw] [--obs-only]
set -euo pipefail
cd /exports/para-lipg-hpc/mdmanurung/cellprograms/benchmarks/biomni_replication
exec /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-py/bin/python -u 01_pseudobulk.py --cohort "$@"
