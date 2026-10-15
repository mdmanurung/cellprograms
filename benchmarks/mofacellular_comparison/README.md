# MOFAcellulaR comparison: scripts and reproduction

Scripts that produce the real-data results in
`vignettes/articles/mofacellular_comparison.Rmd` (the Stephenson COVID-19
PBMC comparison of `cellprograms` against MOFAcellulaR, plus the
artifact-exclusion sensitivity rerun). **No data are committed** — the
pseudobulk, fitted models, and caches are large and are rebuilt from public
data.

## Data

The Stephenson et al. (2021) COVID-19 PBMC cohort, processed to a donor x
gene pseudobulk matrix per cell type (12 cell types, 72-117 donors per
type). Source: ArrayExpress **E-MTAB-8731**. The scripts expect a pseudobulk
object at `data/stephenson_pb/semi_base.rds` — a list with `cts` (per
cell-type `list(Y = donor x gene count matrix, ...)`) and `donor_meta`
(donor-level metadata with `donor_id`, `site`, `sex`,
`Status_on_day_collection_summary`). See
`benchmarks/cmp_detection/scripts/real_data.R` for how this object is
consumed and `vignettes/articles/real_data_case_study.Rmd` for the cohort
description.

## Run order

Set `PB_DIR` (directory holding `stephenson_pb/semi_base.rds`) and
`OUT_DIR` (where fits/caches go; default `results/`). Requires the
`cellprograms` package plus `MOFAcellulaR` and `MOFA2` (with a working
`mofapy2` via reticulate) for steps 3 and 5.

```bash
# 1. Default masked fit (274-gene vdj + genotype_marker blacklist), ~25 min
PB_DIR=data OUT_DIR=results Rscript scripts/01_fit_masked.R

# 2. Cross-cell-type coordination on the default fit, ~1 min
Rscript scripts/02_cmp.R results/stephenson_K30_masked_fit.rds \
  results/stephenson_K30_cmp.rds

# 3. MOFAcellulaR under parity (same HVG panels + blacklist veto), ~5 min
Rscript scripts/03_mofa.R results/stephenson_K30_masked_fit.rds \
  results/stephenson_mofa.hdf5

# 4. Sensitivity refit with exclude_artifacts = TRUE (313-gene blacklist), ~25 min
Rscript scripts/04_fit_excl.R results/stephenson_K30_masked_fit.rds \
  results/stephenson_K30_masked_fit_excl.rds

# 5. CMP on the sensitivity fit, ~1 min
Rscript scripts/02_cmp.R results/stephenson_K30_masked_fit_excl.rds \
  results/stephenson_K30_cmp_excl.rds

# 6. MOFAcellulaR on the sensitivity fit, ~5 min
Rscript scripts/03_mofa.R results/stephenson_K30_masked_fit_excl.rds \
  results/stephenson_mofa_excl.hdf5

# 7. Build the two vignette caches (VIEW_DIMS from the step-3/6 logs)
Rscript scripts/05_build_cache.R results/stephenson_K30_masked_fit.rds \
  results/stephenson_K30_cmp.rds results/stephenson_mofa.hdf5 \
  ../../vignettes/articles/mofacellular_stephenson_cache.rds
VIEW_DIMS=846,433,286,621,376,515,534,250,595,419,874,822 \
  Rscript scripts/05_build_cache.R results/stephenson_K30_masked_fit_excl.rds \
  results/stephenson_K30_cmp_excl.rds results/stephenson_mofa_excl.hdf5 \
  ../../vignettes/articles/mofacellular_stephenson_cache_excl.rds
```

The vignette knits with or without the caches: a guard chunk skips the
real-data results and renders the run-code chunk when a cache is absent, so
the article always builds.

**Exact provenance:** `scripts/as-run/` preserves the three scripts exactly as
executed for the sensitivity arm (including the original machine paths); the
parameterized scripts above are the maintained equivalents.

## Outputs

| File | Content |
|---|---|
| `stephenson_K30_masked_fit.rds` | default masked fit (284 programs) |
| `stephenson_K30_cmp.rds` | CMP result (24 CMPs, 598 edges) |
| `stephenson_mofa.hdf5` | trained MOFA model, default setting |
| `stephenson_K30_masked_fit_excl.rds` | sensitivity fit (278 programs) |
| `stephenson_K30_cmp_excl.rds` | sensitivity CMP (25 CMPs, 516 edges) |
| `stephenson_mofa_excl.hdf5` | trained MOFA model, sensitivity setting |
| `mofacellular_stephenson_cache.rds` | vignette cache, default setting |
| `mofacellular_stephenson_cache_excl.rds` | vignette cache, sensitivity setting |

The two caches are committed under `vignettes/articles/` so the article
reproduces exactly without rerunning the fits.

## Notes

- The simulated-data half of the vignette is self-contained
  (`simulate_pseudobulk(scenario = "mixed")`) and needs none of these
  scripts.
- `VIEW_DIMS` is the per-view gene count after MOFAcellulaR's QC chain,
  printed by step 3/6; it is recorded in the cache for the vignette's QC
  note and does not affect any computation.
- MOFAcellulaR's QC flags complexity-artifact factors (4/8 default, 5/7
  sensitivity); these are a normalization/coverage property of the data,
  reported honestly in the vignette rather than removed.
