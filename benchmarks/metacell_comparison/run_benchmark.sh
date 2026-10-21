#!/bin/bash
# Orchestrator: runs all benchmark stages with reticulate-safe session separation.
# Usage: bash run_benchmark.sh
set -uo pipefail
export R_LIBS=/workspace/.Rlib
BENCH=/workspace/repo/benchmarks/metacell_comparison
WORK=/workspace/bench_mc
FX=$WORK/fixtures
MEM=$WORK/memberships
RES=$BENCH/results
LOGS=$WORK/logs
mkdir -p "$FX" "$MEM" "$RES" "$LOGS"
SCALES="5k 15k"
SEEDS="1 2 3"

echo "=== Stage 1: fixtures ==="
for s in $SCALES; do for seed in $SEEDS; do
  Rscript $BENCH/benchmark_fixture.R "$FX" "$s" "$seed" > "$LOGS/fixture_${s}_${seed}.log" 2>&1 || { echo "FIXTURE FAIL $s $seed"; exit 1; }
done; done
tail -n +1 "$LOGS"/fixture_*.log

echo "=== Stage 2a: SuperCell memberships (pure R) ==="
for s in $SCALES; do for seed in $SEEDS; do
  Rscript $BENCH/benchmark_memberships_supercell.R "$FX" "$MEM" "$s" "$seed" > "$LOGS/mem_sc_${s}_${seed}.log" 2>&1 || echo "SC MEM FAIL $s $seed (continuing)"
done; done
grep -h "SUPERCELL" "$LOGS"/mem_sc_*.log

echo "=== Stage 2b: SEACells memberships (basilisk, own sessions) ==="
for s in $SCALES; do for seed in $SEEDS; do
  Rscript $BENCH/benchmark_memberships_seacells.R "$FX" "$MEM" "$s" "$seed" > "$LOGS/mem_sea_${s}_${seed}.log" 2>&1 || echo "SEA MEM FAIL $s $seed (continuing)"
done; done
grep -h "SEACELLS" "$LOGS"/mem_sea_*.log

echo "=== Stage 3: metrics ==="
for m in supercell seacells; do for s in $SCALES; do for seed in $SEEDS; do for N in 50 100 200 400; do
  if [ -f "$MEM/membership_${m}_${s}_s${seed}_N${N}.csv" ]; then
    Rscript $BENCH/benchmark_metrics.R "$FX" "$MEM" "$s" "$seed" "$m" "$N" "$RES/metrics_${m}_${s}_s${seed}_N${N}.csv" > "$LOGS/metrics_${m}_${s}_${seed}_${N}.log" 2>&1 || echo "METRICS FAIL $m $s $seed $N"
  fi
done; done; done; done; done
grep -h "METRICS" "$LOGS"/metrics_*.log | head -60

echo "=== Stage 4: recovery ==="
for m in supercell seacells; do for s in $SCALES; do for seed in $SEEDS; do for N in 50 100 200 400; do
  if [ -f "$MEM/membership_${m}_${s}_s${seed}_N${N}.csv" ]; then
    Rscript $BENCH/benchmark_recovery.R "$FX" "$MEM" "$s" "$seed" "$m" "$N" "$RES/recovery_${m}_${s}_s${seed}_N${N}.csv" > "$LOGS/recovery_${m}_${s}_${seed}_${N}.log" 2>&1 || echo "RECOVERY FAIL $m $s $seed $N"
  fi
done; done; done; done; done
grep -h "RECOVERY" "$LOGS"/recovery_*.log | head -60

echo "=== Stage 5: report ==="
Rscript $BENCH/benchmark_report.R "$RES" > "$LOGS/report.log" 2>&1
RC=$?
cat "$LOGS/report.log"
echo "BENCH EXIT: $RC"
exit $RC
