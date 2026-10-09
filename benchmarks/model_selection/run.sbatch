#!/bin/bash
#SBATCH --job-name=ms_bench
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/ms.%A_%a.log
# usage (from benchmarks/model_selection): TASKS=tasks_param.txt OFFSET=0 sbatch --array=1-100%2 run.sbatch
# MaxArraySize is 125 here, so a 450-line task file is submitted in chunks: OFFSET=0,100,200,300,400 (task line = OFFSET + array index).
# One array task = one dataset, all arms (paired). Optional env: AMP (pilot override), OUTDIR, ARMS.
HERE=${SLURM_SUBMIT_DIR:-$(pwd)}
ENVBIN=${ENVBIN:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin}
export R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p results/logs
read BASE SCEN REP < <(sed -n "$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))p" "${TASKS:?set TASKS}")
"$ENVBIN/Rscript" "$HERE/run_rep.R" $BASE $SCEN $REP "${AMP:-}" "${OUTDIR:-}" "${ARMS:-}"
