#!/bin/bash
#SBATCH --job-name=real_pa_pmd
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/real_pa_pmd_%j.log
# Usage: cd benchmarks/biomni_replication && sbatch run_real.sbatch
# (#SBATCH --output is relative to the submit dir). Overridable: ENVBIN.
HERE=${SLURM_SUBMIT_DIR:-$(pwd)}
ENVBIN=${ENVBIN:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin}
export R_LIBS_USER=/nonexistent OMP_NUM_THREADS=2
"$ENVBIN/Rscript" "$HERE/06_real_pa_vs_pmd.R"
