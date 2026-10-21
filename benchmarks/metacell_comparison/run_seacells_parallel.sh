#!/bin/bash
# Parallel SEACells driver: one R session per seed (each owns its basilisk Python;
# separate processes, no reticulate conflict). Runs all 4 N targets per seed.
# Usage: bash run_seacells_parallel.sh <scale> <seed1> [seed2] ...
set -u
export R_LIBS=/workspace/.Rlib
BENCH=/workspace/repo/benchmarks/metacell_comparison
FX=/workspace/bench_mc/fixtures
MEM=/workspace/bench_mc/memberships
LOGS=/workspace/bench_mc/logs
mkdir -p "$MEM" "$LOGS"

SCALE=$1; shift

for seed in "$@"; do
  (
    lock="$LOGS/.lock_sea_${SCALE}_${seed}"
    if mkdir "$lock" 2>/dev/null; then
      Rscript $BENCH/benchmark_memberships_seacells.R "$FX" "$MEM" "$SCALE" "$seed" \
        > "$LOGS/mem_sea_${SCALE}_${seed}.log" 2>&1
      rc=$?
      rmdir "$lock"
      echo "SEED DONE $SCALE s$seed rc=$rc"
    else
      while [ -d "$lock" ]; do sleep 20; done
      echo "SEED WAITED $SCALE s$seed"
    fi
  ) &
done
wait
echo "SEACELLS PARALLEL DONE $SCALE"
