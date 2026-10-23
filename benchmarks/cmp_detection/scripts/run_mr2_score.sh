#!/bin/bash
## run_mr2_score.sh — MR-2 scoring driver (K30 vs cNMF, paired).
## Runs the frozen caller (n_agg_noanc via run_cand_v2.R) on each arm's cached
## fits, then the paired bootstrap (v2_score.R). Scores whatever datasets exist
## on the common paired subset; n_datasets in the output discloses the count.
## Usage: bash run_mr2_score.sh <fitroot> <outroot> <cores> [seeds]
##   fitroot: dir containing <arm>/semi/<scenario>/sNNNN.rds
##   outroot: dir for evals/ output (CMP_ROOT of the checkout that owns fitroot)
set -euo pipefail
FITROOT="$1"; OUTROOT="$2"; CORES="$3"; SEEDS="${4:-1051:1075}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CMP_ROOT="$OUTROOT"
export CMP_FITROOT="$FITROOT"
export CMP_SEEDS="$SEEDS"
export CMP_TAG="_mr2"
export R_LIBS="${R_LIBS:-/workspace/.Rlib}"

cd "$HERE"
for ARM in K30 cnmf; do
  N=$(find "$FITROOT/$ARM" -name "s*.rds" 2>/dev/null | wc -l)
  if [ "$N" -eq 0 ]; then echo "skip arm $ARM (no fits)"; continue; fi
  echo "== arm $ARM: $N fits =="
  Rscript run_cand_v2.R n_agg_noanc "$ARM" heldout2 "$CORES" semi
done

## paired scoring on the common subset (K30 = reference arm)
if [ -f "$OUTROOT/evals/heldout2_mr2/K30/n_agg_noanc.rds" ]; then
  REF="$OUTROOT/evals/heldout2_mr2/K30/n_agg_noanc.rds"
  ARGS="$REF=K30"
  for CN in cnmf; do
    F="$OUTROOT/evals/heldout2_mr2/$CN/n_agg_noanc.rds"
    [ -f "$F" ] && ARGS="$ARGS $F=$CN"
  done
  Rscript v2_score.R "$OUTROOT/evals/heldout2_mr2/mr2_paired.csv" "$SEEDS" "$ARGS"
  echo "paired scoring written"
else
  echo "no K30 eval rds; skip paired scoring"
fi
