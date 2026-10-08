#!/bin/bash
#SBATCH --job-name=fit_cp
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/fit_cp.%j.log
# usage (from benchmarks/biomni_replication; #SBATCH --output is relative to the submit dir):
#   sbatch 03_fit_cp.sbatch <cohort> <K> <biomni|repo> <out_rds>
# (medium partition caps at 4 h; L0 fit took ~20 min in the reference -- use --partition=all --time=.. for longer)
# env lives in the archive conda envs_dir (see ~/.condarc); override with ENVBIN=...
HERE=${SLURM_SUBMIT_DIR:-$(pwd)}
ENVBIN=${ENVBIN:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin}
export R_LIBS_USER=/nonexistent
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
"$ENVBIN/Rscript" "$HERE/03_fit_cp.R" "$@"
