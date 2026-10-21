#!/bin/bash
# Run the full program-recovery grid for one scale on this machine.
# Usage: run_recovery_all.sh <scale> <seeds_csv> <fx_dir> <mem_dir> <out_dir> <log_dir> [methods_csv]
set -u
SCALE=$1; SEEDS=$2; FX=$3; MEM=$4; OUT=$5; LOGS=$6
METHODS=${7:-supercell,seacells}
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT" "$LOGS"
for SEED in $(echo "$SEEDS" | tr ',' ' '); do
  for METHOD in $(echo "$METHODS" | tr ',' ' '); do
    for N in 50 100 200 400; do
      CSV="$OUT/recovery_${METHOD}_${SCALE}_s${SEED}_N${N}.csv"
      if [ -s "$CSV" ]; then echo "SKIP $CSV"; continue; fi
      echo "RUN recovery $METHOD $SCALE s$SEED N$N"
      Rscript "$HERE/benchmark_recovery.R" "$FX" "$MEM" "$SCALE" "$SEED" "$METHOD" "$N" "$CSV" \
        > "$LOGS/recovery_${METHOD}_${SCALE}_s${SEED}_N${N}.log" 2>&1 \
        || echo "FAIL recovery $METHOD $SCALE s$SEED N$N"
    done
  done
done
echo "RECOVERY GRID DONE scale=$SCALE"
