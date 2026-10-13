# cellprograms 0.2.0

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
