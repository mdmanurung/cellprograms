#!/usr/bin/env bash
# 01_download_data.sh -- download DIVAS COVID-19 case-study inputs + precomputed results
#
# Data sources:
#   - Input blocks + sample metadata: ByronSyun/DIVAS_COVID19_CaseStudy (GitHub, main branch)
#   - Precomputed DIVAS 6-block results: Zenodo record 17430294
#     (Su et al. COVID-19 multi-omics case study; results from run_combined_6block_divas.R,
#      DIVASmain(nsim=500, colCent=FALSE, rowCent=TRUE, seed=123))
#
# Usage: DATA_DIR=<dir> bash 01_download_data.sh
#   DATA_DIR defaults to ./data ; set to a shared/checkpoint location for large runs.
set -euo pipefail

DATA_DIR="${DATA_DIR:-data}"
INPUT_DIR="$DATA_DIR/input"
mkdir -p "$INPUT_DIR"

GH_BASE="https://raw.githubusercontent.com/ByronSyun/DIVAS_COVID19_CaseStudy/main/scRNA_celltyist_analysis/DIVAS_run/divas_input_T1T2combined"
FILES=(Combined_CD4_T_combined.csv Combined_CD8_T_combined.csv Combined_CD14_Monocyte.csv
       Combined_NK.csv Combined_proteomics.csv Combined_metabolomics.csv
       final_sample_info.csv data_summary.csv)

for f in "${FILES[@]}"; do
  echo "== $f"
  curl -sfL "$GH_BASE/$f" -o "$INPUT_DIR/$f"
done

ZENODO="https://zenodo.org/api/records/17430294/files"
echo "== divas_results_combined_6block_renamed.rds (419 MB, this is the slow one)"
curl -sfL "$ZENODO/divas_results_combined_6block_renamed.rds/content" \
     -o "$DATA_DIR/divas_results_combined_6block_renamed.rds"

echo "== downloaded files:"
ls -la "$INPUT_DIR" "$DATA_DIR"/*.rds
echo "== expected dims (data_summary.csv):"
cat "$INPUT_DIR/data_summary.csv"
echo "DONE"
