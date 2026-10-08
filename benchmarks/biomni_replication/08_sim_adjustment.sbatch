#!/bin/bash
#SBATCH --job-name=sim_adj
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-20
#SBATCH --output=results/logs/sim_adj.%A_%a.log
# usage (from benchmarks/biomni_replication): sbatch 08_sim_adjustment.sbatch   (one rep per array task, 4 scenarios x 2 arms ~1.5 min each)
HERE=${SLURM_SUBMIT_DIR:-$(pwd)}
ENVBIN=${ENVBIN:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin}
export R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p results/sim_adj3
"$ENVBIN/Rscript" "$HERE/08_sim_adjustment.R" $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID results/sim_adj3/rep${SLURM_ARRAY_TASK_ID}.csv
