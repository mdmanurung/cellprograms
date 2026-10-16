# cellprograms 0.2.0.9000 (development)

## Audit remediations

* `find_coordinated_programs()`: missing donor scores (NA) are no longer
  silently converted into top ranks. Rank transformation and strata
  residualization keep NA donors as NA, cross-cell-type correlations use
  pairwise-complete donors, edge `n` reports the pairwise-complete donor
  count for each program pair, and a warning lists programs with missing
  scores. Previously, `rank()`'s default `na.last=TRUE` manufactured fake
  top-scoring donors, which either destroyed detection power (donor-specific
  missingness) or inflated significance (missingness shared across programs);
  a planted pair with 20% missingness could be missed entirely with no
  warning. Regression tests: `tests/testthat/test-cmp-missingness.R`.

# cellprograms 0.2.0

## Stage B: frozen preprocessing and cell-level projection

* `freeze_preprocessing()` (MR-0): resolve the feature panel and adjustment
  decisions once into a content-addressed specification
  (`cell_program_preproc_spec`, sha256 digest); refit any cohort carrying the
  same genes with `fit_celltype_programs(preprocessing_spec = ...)`, which
  replays the spec fail-closed (conflicting explicit arguments are refused,
  and the spec's digest is recorded in the fit's preprocessing provenance).
* `project_programs_to_cells()` (MR-1): project fitted programs onto
  individual cells. Query expression is adjusted with the fit's recorded
  preprocessing (residualization coefficients, centering, scale), activities
  are inferred against the frozen gene loadings via
  `flashier.utils::project_factor_activity()`, and raw activities are
  canonicalized with the same RMS/sign multipliers as the donor scores.
  Requires a canonicalized fit with declared `feature_scale`; fits with
  residualized covariates require `donors` + `sample_metadata`, and the
  rebuilt covariate design must match the fit's design columns exactly.
* `fit_celltype_programs()` gains `feature_scale` (declared measurement
  units, recorded in the fit manifest and required for projection) and
  `preprocessing_spec` (MR-0 replay); fits now record a `preproc_replay`
  block per cell type so the exact adjustment can be replayed on query data.
* New roadmap: `reports/STAGE_B_ROADMAP.md` (MR-0..MR-4 spec, staged main-run
  costs, registration drafts, gates).

## Backend synchronization (Stage A)

* Ordinary EBMF fits now delegate to `flashier.utils::fit_ebmf()` and
  canonicalization uses `flashier.utils::standardize_factors(scaling = "rms")`,
  with numerical parity verified to machine precision against the previous
  direct `flashier::flash()` backend. The public workflow, output fields, and
  scientific results are unchanged. `flashier`, `ebnm`, and `flashier.utils`
  move to `Imports` (R >= 4.3 required); the fixed-covariate path
  (`covariate_mode = "fixed"`) and degenerate fits remain on the legacy
  backend. Set `options(cellprograms.use_flashier_utils = FALSE)` to force the
  legacy backend. See `reports/SYNC_ROADMAP.md` for the version pair, parity
  evidence, and rollback plan.

## New features

* `find_coordinated_programs()` and the CMP toolkit: coordinated multicellular
  program detection from donor-level program scores with within-strata
  permutation FWER control, plus `cmp_scores()`, `ambient_program_scores()`,
  `annotate_cmp_lr()`, `annotate_cmp_themes()`, `cmp_gene_overlap()`,
  `cmp_context_association()`, `annotate_cmp_tf()` (CollecTRI footprints via
  `prepare_collectri()`), `cmp_triage()`, and `cmp_nuisance_genes()`.
* `program_sharing()`: post-hoc shared-vs-specific assessment of matched
  programs across cell types with 'mashr' — permutation-calibrated loading
  matching, greedy agglomeration into disjoint matched groups, CorShrink-style
  null-correlation estimation (`z_thresh` configurable), and per-gene lfsr
  calls. Fits now store loading posterior SDs (`F_psd`/`L_psd`).
* Artifact-gene masking: `artifact_gene_sets` (V(D)J, genotype-marker, sex
  sets), `artifact_genes()`, `flag_artifact_genes()`, and fit-time
  `exclude_genes` in `fit_celltype_programs()` (applied before variable-gene
  selection).

## Bug fixes and robustness

* `program_sharing()` handles degenerate real-data fits: numerically
  degenerate programs (noise-floor loadings; fixed factors with collapsed
  posterior SDs) are excluded from matching, and zero-SD point-laplace
  posteriors are floored relative to each condition's loading scale instead of
  tripping mashr's zero check.
* `simulate_pseudobulk()` no longer fails with custom `cell_types` under the
  "mixed" and "same_genes_independent" scenarios (shared-program objects are
  created by the first present member instead of a hardcoded cell type).

## Documentation and infrastructure

* Full roxygen2 documentation: `man/` pages for all exported functions with
  runnable examples; NAMESPACE is now roxygen-managed.
* R-CMD-check CI (ubuntu + macOS) and a pkgdown site skeleton.
* Quickstart vignette and a precomputed real-data case-study article.

## Benchmarks

* `benchmarks/cmp_detection/`: pre-registered v1.1/v1.2 confirmations of the
  frozen CMP caller (`max_factors = 30` upstream; held-out seeds 1051-1100,
  all 7 criteria passed) and dev verdicts for four challenger arms
  (joint-subspace rejected; two-stage identical to winner; Cauchy closest
  miss; EB shrinkage no-op).
