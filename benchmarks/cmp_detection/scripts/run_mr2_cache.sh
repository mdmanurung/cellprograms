#!/bin/bash
## MR-2 fit-cache driver (PREREG_mr2-cnmf-v1.md): one arm, one seed block.
## Usage: run_mr2_cache.sh <arm: K30|cnmf> <seed_from> <seed_to> <cores> [outroot]
## Env: CMP_ROOT (default: repo checkout two levels up), SEMI_DIR, CNMF_PY.
## Resumable: cache_one skips existing non-empty s%04d.rds files.
set -u
export PATH="/opt/conda/bin:$PATH"
ARM=${1:?arm}; S0=${2:?seed_from}; S1=${3:?seed_to}; CORES=${4:?cores}
OUTROOT=${5:-}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
export CMP_ROOT="${CMP_ROOT:-$REPO}"
export SEMI_DIR="${SEMI_DIR:-$HERE/../data/stephenson_pb}"
export CNMF_PY="${CNMF_PY:-python3}"
SCENS="null,cmp_distinct,cmp_mixed,genes_only,context,ambient,lr,weak"
R_LIBS="${R_LIBS:-$REPO/../../.Rlib}"
Rscript "$HERE/fit_cache_v2.R" "$ARM" semi "$SCENS" "$S0" "$S1" 0.875 "$CORES" $OUTROOT
