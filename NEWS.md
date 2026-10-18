# cellprograms 0.2.0.9000 (development)

## New features (2026-10-18 feature-completion segment)

* `fit_celltype_programs()` gains `backend = c("flashier", "cebmf_torch")`.
  The `cebmf_torch` backend runs the pure-PyTorch cEBMF implementation
  (pinned upstream `william-denault/cebmf_torch@34abac5`) through the
  `cebmf_ext` adapter, with a registered numerical guard set (sebetahat floor
  1e-4, tau precision cap 1e8, explicit seeding, seed-family fallback with
  diagnostics) and a strict mapping contract: the bridge never manufactures
  flashier-only quantities, so `elbo`, `L_psd`, `F_psd`, and `F_lfsr` are NULL
  on this backend and `program_summary()` reports `elbo = NA`. Fixed-mode
  covariates are refused on this backend (no fixed-loadings path); use the
  default `covariate_mode = "residualize"`. Python resolution:
  `options(cellprograms.cebmf_python_env = ...)` (python binary or venv with
  `cebmf_ext` importable), else basilisk (pinned `torch==2.7.1` CPU — the pin
  is load-bearing, 2.14.1+cpu diverges on inputs where 2.7.1+cpu is green),
  else a fail-closed error. No flashier-style fallback ladder: backend
  failures abort loudly per cell type.
* `fit_celltype_programs()` gains `var_floor = "none" | "vignette" | <positive
  scalar>`: fix the residual SD instead of estimating it. `"vignette"` uses
  the flashier single-cell vignette's floor, `S = sd(log(Poisson(1/n) + 1))`
  with `n` the cell type's donor count, drawn once per cell type at a fixed
  local seed (deterministic, fit-seed-independent). A `var_floor`-derived S is
  not mappable through `flashier.utils::fit_ebmf()`, so on the fast path the
  fit reroutes through `flashier::flash()` with a warning; the cebmf_torch
  backend passes S through as KNOWN noise. The default `"none"` preserves
  existing behavior exactly. Calibration evidence: CAND-13 demonstrator
  (`reports/audit/E_programs/`).
* Covariate-audit patches (R2/F1-F4): gene-drop blocks now subset the
  projection-replay `center`/`scale` vectors to retained genes (F1);
  `covariate_mode = "fixed"` with `center = FALSE` is refused even when the
  combination arrives through a frozen preprocessing spec (F3); covariate
  terms referencing within-cell-type-constant variables — including
  interaction terms — are dropped instead of crashing `model.matrix()` (F4);
  aliased rank-deficient replay designs zero their NA coefficients at
  projection instead of propagating NaN activities (F2).

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

## Remediation queue (2026-10-18): CAND-01, CAND-05, CAND-07, D-2a, CAND-03+D-3b, CAND-12

* `assess_program_stability()` (CAND-01, P1): subsample refit failures are now
  counted and reported (`fit$stability$n_failed`/`n_success`, warning on any
  failure); `recovery_freq` is defined over SUCCESSFUL refits — previously the
  silent denominator was `n_boot`, so a 30% refit-failure rate read as a 30%
  recovery loss with zero user-facing signal, and an all-fail run was
  indistinguishable from zero stability. All-fail now yields NA with a
  distinct warning; per-program `n_hit` added.
* `find_coordinated_programs()` (CAND-05, P2): CMP member orientation now
  propagates along the linking edges that formed the CMP (product of edge
  signs along the path from the reference member). Chain-joined members
  previously inherited the sign of the direct reference–member edge from the
  full edge table — which can be a non-linking noise edge, i.e. a random +/-
  that propagated into `cmp_scores()`, `annotate_cmp_tf()` and `cmp_triage()`
  direction calls.
* `program_sharing()` (D-2a/D-2c, P1): `matches=` references are validated
  (unknown programs abort informatively instead of crashing in the mashr
  layer) and every forced-matches call warns that the "shared" label counts
  top genes ACTIVE (lfsr < thresh) in >= 2 conditions — which independent
  loadings can also reach at high loading SNR (audit D-2a: 30/30 false
  "shared" for forced independent pairs at loading SE 0.125, while mashr's
  own mixture put 96% weight on the identity pattern, p_all_shared ~ 0.001);
  `p_all_shared` is the diagnostic. The classification rule itself is
  unchanged pending a registered benchmark of any replacement.
* Matching-step documentation (CAND-07, P1): the gene-permutation null's
  exchangeability caveat is now stated where matching q-values are produced
  (mirroring `cmp-interpret.R`): block-correlated genes inflate false matches
  (audit D-1: 67.0% type-I under block equicorrelation 0.55 vs 2.0% iid). The
  null itself is UNCHANGED — a registered mini-benchmark (null calibration at
  block equicorrelation + planted-truth power) precedes any null change.
* `as_cell_program_data()` (CAND-03, P2): rejects ALL non-finite entries
  (NA/NaN/Inf) with a count-bearing abort; Inf previously passed the NA-only
  check and crashed validate/fit with an opaque
  "missing value where TRUE/FALSE needed".
* `.greedy_match()` (D-3b, P2): non-finite similarities are demoted to
  "unmatched" with a warning instead of crashing on `max()` of an NA-bearing
  matrix.
* `fit_celltype_programs()` + `assess_program_stability()` (CAND-12, P2): the
  full frozen `preprocessing_spec` is stored on the fit (not just its digest)
  and stability refits REPLAY it through the MR-0 fail-closed checks;
  `feature_scale` joined the MR-0 conflict list and is passed through in
  legacy-path refits (previously silently "unknown" in refit manifests).

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
