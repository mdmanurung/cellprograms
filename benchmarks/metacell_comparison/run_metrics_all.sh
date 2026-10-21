#!/bin/bash
# Run the full metrics grid for one scale on this machine.
# Usage: run_metrics_all.sh <scale> <seeds_csv> <fx_dir> <mem_dir> <out_dir> <log_dir> [methods_csv]
# Example: run_metrics_all.sh 5k 1,2,3 /workspace/bench_mc/fixtures /workspace/bench_mc/memberships /workspace/bench_mc/results /workspace/bench_mc/logs
set -u
SCALE=$1; SEEDS=$2; FX=$3; MEM=$4; OUT=$5; LOGS=$6
METHODS=${7:-supercell,seacells}
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT" "$LOGS"
for SEED in $(echo "$SEEDS" | tr ',' ' '); do
  for METHOD in $(echo "$METHODS" | tr ',' ' '); do
    for N in 50 100 200 400; do
      CSV="$OUT/metrics_${METHOD}_${SCALE}_s${SEED}_N${N}.csv"
      if [ -s "$CSV" ]; then echo "SKIP $CSV"; continue; fi
      echo "RUN metrics $METHOD $SCALE s$SEED N$N"
      Rscript "$HERE/benchmark_metrics.R" "$FX" "$MEM" "$SCALE" "$SEED" "$METHOD" "$N" "$CSV" \
        > "$LOGS/metrics_${METHOD}_${SCALE}_s${SEED}_N${N}.log" 2>&1 \
        || echo "FAIL metrics $METHOD $SCALE s$SEED N$N"
    done
  done
done
echo "METRICS GRID DONE scale=$SCALE"
