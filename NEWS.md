# cellprograms 0.2.0.9000 (development)

## Bug fixes and hardening (2026-10-19 implementation-review segment)

Full implementation review follow-through: 2 critical correctness fixes, 12
major fixes, and 50 minor fixes. Regression coverage for every fixed finding
lives in `tests/testthat/test-review-fixes.R` and
`tests/testthat/test-artifact-genes.R`.

### Critical

* `GetMetacellSEACells()`: multi-group soft-assignment matrices are now
  combined block-diagonally (one block per group) instead of row-binding
  group-local column indices, which crashed or silently misassigned cells to
  metacells when more than one group was processed.
* `cmp_context_association()`: the context term is selected exactly through
  the model-matrix term assignment instead of a `^x` prefix grep that also
  captured adjustment covariates whose names start with "x" (e.g. `xsite`).
  Covariates named `y` or `.ctx` are refused (reserved internal names).

### Major

* `fit_celltype_programs()`: within-cell-type constant genes no longer crash
  the gene filter under `scale = TRUE`; their non-finite column mass is
  treated as degenerate and the genes are dropped.
* `fit_celltype_programs(backend = "cebmf_torch")` warns when a cell type's
  fit exhausts `maxit` without converging (the flag was stored but never
  surfaced).
* `assess_program_stability()` warns and skips cell types with fewer than 5
  donors instead of crashing during resampling.
* `plot_program_scores()` handles per-cell-type observation sets (union of
  donors, NA where a cell type was not observed) instead of requiring equal
  row counts; z-scoring is NA-aware.
* `plot_sharing_spectrum()` rewritten: significance class (shared / n.s. /
  underpowered) drives a fixed colorblind-safe palette, q-values are carried
  in the pair labels, and >9 pairs no longer break the color scale.
* `plot_artifact_flags()` detects any `*_frac` column set (not just the
  built-in screens) and accepts CMP-shaped flag tables.
* `union_panel_loadings()` leaves zero-overlap or constant panel genes NA
  (with warnings) instead of aborting, and replays rank-deficient fixed
  covariate designs via least squares instead of crashing in `solve()`.
* `project_programs_to_cells()` refuses a `feature_scale` declaration that
  contradicts the fit's recorded units (previously silently accepted).
* `flag_artifact_genes()`: the V(D)J screen now also covers J/D segments,
  IGLL-family genes, and JCHAIN.
* Metacells: `max_iter` documentation corrected (wrapped SEACells default is
  20); `use_sparse` warns when set (ignored by the wrapped SEACells version);
  per-metacell QC columns are joined into the result metadata.
* Metacell IO coercions use the non-deprecated `CsparseMatrix` form
  throughout: the pinned flashier's `.onLoad` sets
  `Matrix.warnDeprecatedCoerce = 2`, under which the previous `dgCMatrix`
  coercions of `readMM()` output error out (exposed by the first
  `CPMC_INTEGRATION=1` end-to-end run, which now passes).
* `withr` moved from Suggests to Imports: `withr::local_seed()` is called
  unconditionally by the fitting code since the scoped-RNG fix.
* `OptimizeGammaSEACells()` now works under Seurat v5: mcRigor's
  `buildmc` aligns per-metacell sizes through `AggregateExpression`
  colnames, which mangle membership labels (digit-leading values gain a
  "g" prefix; "_" becomes "-") and only undo the digit case — the
  `"mc_<k>"` labels used here never realigned, so every group failed with
  "No cell overlap between new meta data and Seurat object". Membership is
  now integer-coded per gamma column (the same contract as the
  `mcRigor_DETECT` path) and gamma candidates whose membership collapses
  to a single metacell are dropped with a warning (`NA` in `scores`)
  instead of failing the whole group's optimization.

### Minor (selected)

* All internal RNG is now scoped with `withr::local_seed()`; no `set.seed()`
  leaks into the user session.
* `plot_variance_explained()` maps PVEs to programs through the flashbridge
  native-factor mapping and skips (with a warning) cell types whose PVE
  vector cannot be reconciled, instead of positionally truncating.
* `plot_context_association()` omits confidence intervals for associations
  with a zero t statistic (previously infinite whiskers) and labels the
  context term plainly after the `.ctx` internal rename.
* `plot_program_scores()` fails closed when an annotation data.frame lacks
  an `observation_id` column (previously a silent all-NA strip).
* `match_programs()` documents that `space` only affects the subspace angles
  (factor-level matches are always loading-cosine based) and fails closed
  when `space = "scores"` has fewer than 2 shared observations.
* `as_mofacellular_input()` aborts on composite-name collisions under the
  `<cell_type>_<observation_id>` scheme, on duplicate `cell_counts` keys, and
  on zero-observation cell types.
* Metacell IO: missing counts-like layers warn or abort, `dims` are
  validated, and numeric-like cell names abort (they would be mangled by the
  Python bridge).
* CMP machinery: sign-concordance ties and all-NA sign vectors no longer
  score deterministically; `prepare_collectri()` parses integer 0/1 `mor`
  correctly; NA `mor` edges are dropped with a warning;
  `cmp_gene_overlap()` is schema-stable on empty and single-member input.
* Sharing/alignment: degenerate-program NA propagation fixed;
  `sharing_spectrum()` gains an explicit `seed` argument; the variance
  estimate reuses the sharing matrices (single computation, consistent
  panels); principal-angle p-value NA handling no longer biases the
  permutation null.

## New features (2026-10-18 feature-completion segment)

* Metacell construction is folded into `cellprograms` from the retired
  `r_packages/cellprogramsMC` sub-package: `MetacellsByGroupsSEACells()`
  (group-stratified SEACells metacells with mcRigor QC),
  `OptimizeGammaSEACells()` (mcRigor granularity selection), and
  `GetMetacellSEACells()` (result accessor). The basilisk Python environment
  is now created lazily on first use (internal `.seacells_env()`, same
  pattern as the cebmf backend), so `basilisk` stays in Suggests and
  installing `cellprograms` never requires basilisk, Python, Seurat, or
  mcRigor; the exported `seacells_env` object from the sub-package is
  retired (nothing in the project referenced it). `Seurat`, `SeuratObject`,
  `mcRigor`, and `hdWGCNA` are Suggests with explicit fail-closed guards;
  `Matrix` and `methods` move to Imports. Note: the basilisk environment is
  re-registered under `cellprograms`, so the first metacell run after this
  change rebuilds the SEACells conda environment (one-time, needs network).
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

## New features (spec v2 execution segment)

* `assess_cmp_stability()` (CP-2): stability of coordinated multicellular
  programs under SYNCHRONIZED donor resampling. Donors are drawn once per
  replicate (optionally within donor-level `donor_strata`) and the identical
  donor set is applied to every cell type — independent per-type subsampling
  (`assess_program_stability()`) changes donor overlap arbitrarily and cannot
  answer whether the same donors support a CMP. Resampling is subsampling
  without replacement (`subsample_frac`, default 0.8): duplicated donor rows
  would pseudo-replicate the donor-level evidence the permutation null
  evaluates. `donor_col` promotes resampling to the biological donor when it
  differs from `observation_ids` (all of a donor's observations are kept
  together; missing cell types/visits are never synthesized); `visit_col` is
  required when donors contribute > 1 observation, making the
  paired-observation estimand explicit. Each replicate reruns
  `find_coordinated_programs()` with the caller's predeclared `cmp_args` and
  reports per-CMP membership recovery (Jaccard), assignment-ambiguity
  margins, sign concordance up to one global flip, donor-score concordance,
  cell-type survival, and CP-1-style failure accounting (`valid_refits` /
  `failed_refits`, conditional rates, NA on all-fail). When the reference fit
  yields no CMPs, stability is reported as ineligible rather than zero.

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

* `assess_program_stability()` (CAND-01/CP-1, P1): subsample refit failures
  are now counted and reported (`fit$stability$valid_refits`/`failed_refits`,
  warning on any failure) with a per-replicate `replicate_log` (status, error
  class/message, elapsed time, sampled observations and factor counts per
  cell type, matching eligibility). `recovery_freq` keeps its LEGACY meaning
  (unconditional hit fraction over requested iterations, `n_hit / n_boot`);
  the new `conditional_recovery_freq` is defined over SUCCESSFUL refits
  (`n_hit / valid_refits`) — previously the silent denominator was `n_boot`,
  so a 30% refit-failure rate read as a 30% recovery loss with zero
  user-facing signal, and an all-fail run was indistinguishable from zero
  stability. All-fail now yields NA rates with a distinct warning; per-program
  `n_hit` and `n_eligible` added. Subsampled refit objects now carry
  consistent `observation_ids`/`sample_metadata`/`dims` (previously stale
  full-data values).
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

## Cross-cohort replication case study (CP-5, 2026-10-19)

* New evaluated article `vignettes/articles/cross_cohort_replication.Rmd`:
  an end-to-end replication workflow on a three-cohort simulated fixture
  (conserved outcome-linked program, cohort-specific program, nuisance-only
  batch effect) orchestrating existing exports only — no new API. Shows the
  fixed-basis pathway (`project_programs_to_cells()`, learned transform
  applied) versus the de-novo pathway (`freeze_preprocessing()` +
  independent refit + `match_programs()`), with predeclared
  same-sign-as-discovery replication rules, abstention accounting
  (attempted/projected/matched/ineligible/replicated denominators), three
  negative controls (200-fold permuted outcome at null ~0.05; scrambled
  gene IDs abstain loudly; batch-only effect does not replicate as
  biology), and a leave-one-cohort-out check in which the cohort-specific
  program correctly abstains. `confirm_programs_across_cohorts()` is
  deliberately NOT exported: no multi-user need is proven.
