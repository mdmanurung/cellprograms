#!/bin/bash
# usage: make_tasks.sh <param|semi> [rep_offset]  -> tasks_<base>.txt ("base scenario rep" per line)
# Reps per scenario (PREREG): S0 200, S1 50, S2 50, S3 100, S7 50. rep_offset 10000 gives the confirmation seeds.
B=$1; O=${2:-0}; OUT=tasks_${B}${O/#0/}.txt; : > $OUT
for s in "S0 200" "S1 50" "S2 50" "S3 100" "S7 50"; do set -- $s; for r in $(seq 1 $2); do echo "$B $1 $((r+O))" >> $OUT; done; done
echo $OUT $(wc -l < $OUT)
