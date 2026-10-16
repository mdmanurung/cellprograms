# Changelog

## cellprograms 0.2.0

### Stage B: frozen preprocessing and cell-level projection

- [`freeze_preprocessing()`](https://mdmanurung.github.io/cellprograms/reference/freeze_preprocessing.md)
  (MR-0): resolve the feature panel and adjustment decisions once into a
  content-addressed specification (`cell_program_preproc_spec`, sha256
  digest); refit any cohort carrying the same genes with
  `fit_celltype_programs(preprocessing_spec = ...)`, which replays the
  spec fail-closed (conflicting explicit arguments are refused, and the
  spec’s digest is recorded in the fit’s preprocessing provenance).
- [`project_programs_to_cells()`](https://mdmanurung.github.io/cellprograms/reference/project_programs_to_cells.md)
  (MR-1): project fitted programs onto individual cells. Query
  expression is adjusted with the fit’s recorded preprocessing
  (residualization coefficients, centering, scale), activities are
  inferred against the frozen gene loadings via
  [`flashier.utils::project_factor_activity()`](https://mdmanurung.github.io/flashier-utils/reference/project_factor_activity.html),
  and raw activities are canonicalized with the same RMS/sign
  multipliers as the donor scores. Requires a canonicalized fit with
  declared `feature_scale`; fits with residualized covariates require
  `donors` + `sample_metadata`, and the rebuilt covariate design must
  match the fit’s design columns exactly.
- [`fit_celltype_programs()`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)
  gains `feature_scale` (declared measurement units, recorded in the fit
  manifest and required for projection) and `preprocessing_spec` (MR-0
  replay); fits now record a `preproc_replay` block per cell type so the
  exact adjustment can be replayed on query data.
- New roadmap: `reports/STAGE_B_ROADMAP.md` (MR-0..MR-4 spec, staged
  main-run costs, registration drafts, gates).

### Backend synchronization (Stage A)

- Ordinary EBMF fits now delegate to
  [`flashier.utils::fit_ebmf()`](https://mdmanurung.github.io/flashier-utils/reference/fit_ebmf.html)
  and canonicalization uses
  `flashier.utils::standardize_factors(scaling = "rms")`, with numerical
  parity verified to machine precision against the previous direct
  [`flashier::flash()`](https://rdrr.io/pkg/flashier/man/flash.html)
  backend. The public workflow, output fields, and scientific results
  are unchanged. `flashier`, `ebnm`, and `flashier.utils` move to
  `Imports` (R \>= 4.3 required); the fixed-covariate path
  (`covariate_mode = "fixed"`) and degenerate fits remain on the legacy
  backend. Set `options(cellprograms.use_flashier_utils = FALSE)` to
  force the legacy backend. See `reports/SYNC_ROADMAP.md` for the
  version pair, parity evidence, and rollback plan.

### New features

- [`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
  and the CMP toolkit: coordinated multicellular program detection from
  donor-level program scores with within-strata permutation FWER
  control, plus
  [`cmp_scores()`](https://mdmanurung.github.io/cellprograms/reference/cmp_scores.md),
  [`ambient_program_scores()`](https://mdmanurung.github.io/cellprograms/reference/ambient_program_scores.md),
  [`annotate_cmp_lr()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_lr.md),
  [`annotate_cmp_themes()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_themes.md),
  [`cmp_gene_overlap()`](https://mdmanurung.github.io/cellprograms/reference/cmp_gene_overlap.md),
  [`cmp_context_association()`](https://mdmanurung.github.io/cellprograms/reference/cmp_context_association.md),
  [`annotate_cmp_tf()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_tf.md)
  (CollecTRI footprints via
  [`prepare_collectri()`](https://mdmanurung.github.io/cellprograms/reference/prepare_collectri.md)),
  [`cmp_triage()`](https://mdmanurung.github.io/cellprograms/reference/cmp_triage.md),
  and
  [`cmp_nuisance_genes()`](https://mdmanurung.github.io/cellprograms/reference/cmp_nuisance_genes.md).
- [`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md):
  post-hoc shared-vs-specific assessment of matched programs across cell
  types with ‘mashr’ — permutation-calibrated loading matching, greedy
  agglomeration into disjoint matched groups, CorShrink-style
  null-correlation estimation (`z_thresh` configurable), and per-gene
  lfsr calls. Fits now store loading posterior SDs (`F_psd`/`L_psd`).
- Artifact-gene masking: `artifact_gene_sets` (V(D)J, genotype-marker,
  sex sets),
  [`artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/artifact_genes.md),
  [`flag_artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md),
  and fit-time `exclude_genes` in
  [`fit_celltype_programs()`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)
  (applied before variable-gene selection).

### Bug fixes and robustness

- [`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)
  handles degenerate real-data fits: numerically degenerate programs
  (noise-floor loadings; fixed factors with collapsed posterior SDs) are
  excluded from matching, and zero-SD point-laplace posteriors are
  floored relative to each condition’s loading scale instead of tripping
  mashr’s zero check.
- [`simulate_pseudobulk()`](https://mdmanurung.github.io/cellprograms/reference/simulate_pseudobulk.md)
  no longer fails with custom `cell_types` under the “mixed” and
  “same_genes_independent” scenarios (shared-program objects are created
  by the first present member instead of a hardcoded cell type).

### Documentation and infrastructure

- Full roxygen2 documentation: `man/` pages for all exported functions
  with runnable examples; NAMESPACE is now roxygen-managed.
- R-CMD-check CI (ubuntu + macOS) and a pkgdown site skeleton.
- Quickstart vignette and a precomputed real-data case-study article.

### Benchmarks

- `benchmarks/cmp_detection/`: pre-registered v1.1/v1.2 confirmations of
  the frozen CMP caller (`max_factors = 30` upstream; held-out seeds
  1051-1100, all 7 criteria passed) and dev verdicts for four challenger
  arms (joint-subspace rejected; two-stage identical to winner; Cauchy
  closest miss; EB shrinkage no-op).
