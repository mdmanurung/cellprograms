#!/bin/bash
#SBATCH --job-name=adj_share
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/adj_share.%j.log
# usage (from benchmarks/biomni_replication): sbatch 07_adjusted_sharing.sbatch <cohort> <none|cov1,cov2> <out_prefix>
HERE=${SLURM_SUBMIT_DIR:-$(pwd)}
ENVBIN=${ENVBIN:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin}
export R_LIBS_USER=/nonexistent OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
"$ENVBIN/Rscript" "$HERE/07_adjusted_sharing.R" "$@"
