# Experiments index

Top-level index tying each benchmark/study folder to the scripts that ran it
and the figures, reports, and data products it produced. Package source
(`R/`, `tests/`, `man/`) is separate; nothing indexed here runs during
`R CMD check`. Conventions for adding new studies are at the bottom.

## Study index

### `benchmarks/mofacellular_comparison/` — cellprograms vs MOFAcellulaR

- **Question**: does the cellprograms pipeline (masked refit + CMP + mashr
  sharing) recover the same program structure as MOFAcellulaR on real data
  (Stephenson K30), with and without artifact-gene exclusion?
- **Scripts**: `scripts/01_fit_masked.R` ... `scripts/05_build_cache.R`;
  exact executed versions in `scripts/as-run/`.
- **Figures**: `vignettes/articles/figure/cmp-1.png`, `vignettes/articles/figure/ambient-guard-1.png`,
  `vignettes/articles/figure/ambient-guard-2.png`, `vignettes/articles/figure/sharing-1.png`, `vignettes/articles/figure/align-1.png`,
  `vignettes/articles/figure/real-align-1.png`, `vignettes/articles/figure/real-severity-1.png`,
  `vignettes/articles/figure/excl-align-1.png` — all referenced by the rendered article
  `vignettes/articles/mofacellular_comparison.md` (source:
  `mofacellular_comparison.Rmd`).
- **Data products**: committed caches
  `vignettes/articles/mofacellular_stephenson_cache.rds` and
  `vignettes/articles/mofacellular_stephenson_cache_excl.rds` (built by
  `05_build_cache.R`; as-run versions `scripts/as-run/build_cache*.R`).
- **Reports**: narrative lives in the article; sync context in
  `reports/SYNC_ROADMAP.md`.

### `benchmarks/cmp_bench/` — CMP candidate benchmark + sharing/subspace calibration

- **Question**: which candidate program-calling method wins on simulated and
  real cohorts (preregistered in `PREREG.md`); empirical sharing patterns and
  matching calibration on the real K30 fits.
- **Scripts**: harness modules (folder root) + `drivers/` (`run_sharing.R`,
  `run_subspace.R`, `run_b2_b3.R`, smoke/profiling wrappers,
  `summarize_sharing.R`, `breadth_winner.R`).
- **Figures/reports**: most campaign figures and summary CSVs were written to
  the session workspace and were not preserved (see the as-run caveats in
  this folder's README). Committed outputs: `outputs/v_sensitivity.csv` and
  `outputs/matching_calibration.csv` from `run_b2_b3.R` (B2 `z_thresh` sweep,
  B3 matching calibration).
- **Data products**: fit/sharing caches were session-only
  (`/mnt/shared-workspace/cmp/`), not committed.

### `benchmarks/biomni_replication/` — Biomni pseudobulk replication

- **Question**: reproduce Biomni's pseudobulk + program pipeline from Census
  data exactly (validation of our pseudobulk recipe against their execution
  trace).
- **Scripts**: `01_pseudobulk.py` (+ `.sbatch` wrappers), `02*`–`09_*` R
  pipeline, `check_floor.R`, `compare_to_biomni.R`, and validation scripts
  under `data/_validation/`.
- **Figures/reports**: `reports/biomni_cmp_lessons.md`,
  `reports/biomni_transfer/part1_callers_scoring.md`.
- **Data products**: Census-derived pseudobulk caches (`data/`, `mat/`) live
  on the HPC/session store, not committed (recipes and validation results in
  the folder README).

## Package-facing articles fed by experiment work

- `vignettes/articles/real_data_case_study.Rmd` and
  `vignettes/articles/divas_case_study.Rmd`, with committed figures in
  `vignettes/articles/figures/fig_*.png`.
- `vignettes/articles/model_mechanisms.Rmd`.

## Related development reports (not study outputs)

- `reports/final_report.md` (EBMF config-tuning campaign
  `cellprograms-ebmf-config-v1`; its driving scripts predate the
  `benchmarks/` layout), `reports/metaprogram_methods_review.md`,
  `reports/plan_review.md`, `reports/ROADMAP.md`,
  `reports/SYNC_ROADMAP.md`, `reports/plans/`.

## Conventions for new experiments

1. **One study = one folder**: `benchmarks/<study_slug>/`, snake_case slug,
   self-contained (scripts, README, small outputs).
2. **Script naming**:
   - Numbered pipeline scripts in execution order: `01_<step>.R`, `02_<step>.R`
     (cluster wrappers as `01_<step>.sbatch` next to them).
   - Ad-hoc drivers: verb-prefixed snake_case — `run_<thing>.R`,
     `summarize_<thing>.R`, `prof_<thing>.R`, `smoke_<thing>.R`.
   - When a cleaned/numbered script exists, the exact executed version goes in
     `scripts/as-run/` (verbatim copy, path warts included).
3. **README required** in every study folder, with: (a) the question, (b) a
   script table, (c) inputs and where caches come from, (d) outputs and where
   they land, (e) as-run environment caveats (absolute paths, package
   versions, cluster notes).
4. **Outputs**:
   - Commit small figures/tables under the study folder (`outputs/` or
     `figs/`), or as a vignette article when package-facing.
   - Never commit large intermediates (`.rds`/`.hdf5` caches); record their
     path and the script that regenerates them.
   - Mark session-only outputs as such in the README and here.
5. **Load the package, never edit it from experiments**: use
   `pkgload::load_all()` or `library()`; analysis code never goes into `R/`.
6. **Update this index in the same commit** that adds or changes a study
   folder or its outputs.
