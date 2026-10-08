# cellprograms real-data benchmark — organized deliverables

Cross-cohort benchmark of the `cellprograms` R package (per-cell-type empirical-Bayes matrix factorization with the local-first borrowing ladder) against sample-level baselines and comparators (GloScope, PILOT). The package source itself is in the separate `cellprograms/` folder.

## Contents

- `report_benchmark_summary.md` — **start here**: cross-cohort findings, recommendations, caveats.
- `combat/`, `hlca/`, `stephenson/`, `onek1k/` — one folder per cohort, each with:
  - `report_benchmark_<cohort>.md` — per-dataset report (Panels A–E: biology preservation, technical leakage, representation quality, interpretability, practicality).
  - `figures/` — Panel A/B heatmap, R² ladder, block-decomposition mass (COMBAT also: L1 divergence).
  - `tables/` — `eval_all.csv` (all representations × covariates), `program_summary.csv` (per-CT R²), `program_annotations.csv` (top genes + technical screen), `block_structure.csv`, `runtime.csv`, `sharing_spectrum_*.csv`; COMBAT also `l1_divergence_*.csv` and `stability.csv`.
- `cross_cohort/cross_cohort_matches.csv` — COMBAT↔Stephenson program matching (L0 fits, loadings space): 28 matched pairs, 13 strong (|cos| ≥ 0.8).
- `l3_race/` — L3 integration method-selection artifacts (scorecard + state table) that decided block = default backend.

## Cohort configs

| cohort | cells | donors | cell types | ladder run |
|---|---|---|---|---|
| COMBAT (COVID-19 PBMC) | ~1.0M | 50–121/CT | 10 | full (L0/L1-global/L1-groups/L2 + block/pva/evbS) |
| HLCA core (lung atlas) | ~584k | 40–69/CT | 9 | full |
| Stephenson (COVID-19 PBMC) | ~780k | 31–120/CT | 19 | reduced (L0 + L1-groups + block) |
| OneK1K (healthy PBMC) | ~1.3M | 43–300/CT | 9 | reduced (L0 + block) |

Fitted model objects (.rds) are not duplicated here (multi-GB); ask if you want them archived into this folder.
