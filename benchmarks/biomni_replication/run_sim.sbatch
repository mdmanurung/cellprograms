#!/bin/bash
#SBATCH --job-name=sim_pa_pmd
#SBATCH --partition=all
#SBATCH --array=1-20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/sim_%A_%a.log
# Usage (from benchmarks/biomni_replication; #SBATCH --output is relative to the submit dir):
#   sbatch run_sim.sbatch            # one seed per array task; 18 scenario x missing configs each
#   OUT=sim_smoke_v2 sbatch --array=1 run_sim.sbatch --smoke    # smoke test
# Overridable: ENVBIN (conda env bin dir), OUT (results subdir, default sim_v2).
HERE=${SLURM_SUBMIT_DIR:-$(pwd)}
ENVBIN=${ENVBIN:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin}
OUT=${OUT:-sim_v2}
export R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$HERE/results/$OUT" "$HERE/results/logs"
"$ENVBIN/Rscript" "$HERE/05_sim_pa_vs_pmd.R" "$SLURM_ARRAY_TASK_ID" "$HERE/results/$OUT/seed${SLURM_ARRAY_TASK_ID}.csv" "$@"
