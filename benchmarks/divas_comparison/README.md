# DIVAS COVID-19 case study: cellprograms vs DIVAS

Reanalysis of the DIVAS COVID-19 multi-omics case study (Su et al.; data/code:
<https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy>) with cellprograms, and a
head-to-head comparison against the published DIVAS results.

## Design

**Data**: 6 blocks, 114 samples (57 patients x 2 timepoints, T1/T2):
4 scRNA-seq cell-type blocks (CD4+ T, CD8+ T, CD14+ monocytes, NK; common
8,634-gene panel), plasma proteomics (481 features), metabolomics (763 features).
Clinical metadata: WHO ordinal score (wos), severity (mild/moderate/severe),
sex, age.

**cellprograms recipe** (DIVAS-parity, mirrors `run_combined_6block_divas.R`):
`fit_celltype_programs(max_factors = 30, features = "all", center = TRUE,
scale = FALSE, loading_prior = "point_laplace", var_type = 1, backfit = TRUE,
nullcheck = TRUE, seed = 1)` — DIVAS centers features (`rowCent = TRUE`) and
does not scale (`colCent = FALSE`), and runs unadjusted; we match both.
Coordination: `find_coordinated_programs(n_perm = 999, seed = 1)`, marginal
(primary) and within-severity (secondary). Sharing: `program_sharing` on the 4
gene blocks (`n_perm = 199, z_thresh = 1.5`; mashr requires a common feature
space, so proteomics/metabolomics participate in coordination only).

**DIVAS reference**: precomputed results from Zenodo record 17430294
(`divas_results_combined_6block_renamed.rds`; `DIVASmain(nsim = 500, seed = 123)`).

**Comparison**:
1. Structural concordance — |cor| between DIVAS component scores and
   cellprograms program/CMP scores; top-100 loading overlap (Jaccard) for
   matched gene-block programs; taxonomy agreement (DIVAS
   joint/partial/individual vs cellprograms shared/partial/private).
2. Clinical association — component/program association with wos (Spearman)
   and severity (Kruskal), BH-adjusted; cross-validated ridge prediction of
   wos from scores (identical folds for both methods).
3. Unique findings — components/programs with no counterpart (|cor| < 0.3).

**Sensitivities**: T1/T2 samples are paired (same patients); primary analysis
follows the DIVAS combined design, with (a) patient-blocked permutation of
headline CMP calls and (b) a sex+age-residualized refit.

## Layout

- `scripts/01_download_data.sh` — inputs (GitHub) + precomputed DIVAS results (Zenodo)
- `scripts/02_fit_cellprograms.R` — fit + canonicalize all 6 blocks
- `scripts/03_cmp_sharing.R` — coordination + sharing (+ `--sensitivity` refit)
- `scripts/04_extract_divas.R` — extract DIVAS scores/loadings/taxonomy to CSV
- `scripts/05_compare.R` — concordance + clinical comparison + figures
- `results/` — result tables (small CSVs only; large objects live in data/)
- `REPORT.md` — findings

## Reproducing

```bash
DATA_DIR=data bash scripts/01_download_data.sh
Rscript scripts/02_fit_cellprograms.R data results
Rscript scripts/03_cmp_sharing.R data results
Rscript scripts/04_extract_divas.R data results
Rscript scripts/05_compare.R data results
```

Requires R >= 4.2 with cellprograms (plus its flashier/ebnm/mashr deps) and
the DIVAS package only if rerunning DIVAS from scratch (not required: the
precomputed results object is used).
