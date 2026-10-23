# flashier.utils 1.1.0.9000

- Add `activity_threshold`/`activity_direction` to `summarize_factor_activity()`:
  strict greater-than prevalence summaries (`active_n`, `active_proportion`, and
  `active_weighted_proportion` with weights) per factor and group on the resolved
  representation coordinates. Without a threshold the table is unchanged.
- Add `order_within_group = "mds1"` to `plot_factor_activity()` for structure
  plots: display-only within-group ordering by 1D classical MDS of the displayed
  row-sum-normalized proportions (Euclidean distance, first-appearance group
  order, deterministic input-order tie-breaks). Default `"input"` is unchanged.
- Add fresh-native verification oracles for the FND-001 projection manifest
  repair and the FND-002 logical-mask partial fixing, plus a multistart
  `init_seeds` input-order invariance regression; the original delivery
  receipts were session-local and are lost (see HANDOFF).

- Add opt-in `standardize_factors(scaling = "rms")`: rescales each activity
  column to unit root-mean-square and compensates on the effects side,
  preserving the reconstruction exactly. This is a factor-coordinate change
  (units), not a change to the fitted model or its uncertainty; the per-factor
  scale is recorded in the basis manifest as `rms_scale`. Defaults
  (`raw`/`ldf`) are unchanged. Added to support downstream `cellprograms`
  canonicalization parity (SYNC-1A).

- Add native `fit_multistart_ebmf()` with distinct initializer seeds, explicit
  sign constraints, candidate failure/ELBO/digest records, deterministic ties,
  optional candidate retention and callback-required stability replay. Selection
  measures optimization sensitivity and does not establish biological certainty.
- Intercept, initialized and projected views now retain validated per-factor
  support and roles in their basis hashes. Ordinary legacy manifests remain
  unchanged; projection preserves fixed-side support and records estimated priors.
- Zero-factor selection, association, backprojection, network, display and
  stability consumers fail with a clear biological-factor error. Empty extraction
  and summaries remain available, as does default intercept-only backfitting.
- Automatic biological bootstrap replays recorded feature-intercept recipes in
  both orientations; incomplete intercept recipes require an explicit callback.
- Add `fit_initialized_ebmf()` for named sample/feature initialization and
  optional one-side anchoring. Initial IDs are retained, additional IDs avoid
  collisions, and automatic nullchecking only removes new adaptive factors.
  Input hashes are recorded; stability refits require an explicit callback.
- Add opt-in feature intercepts to `fit_ebmf()`. Their sample-side ones are
  protected, offsets remain separate from semantic factors, and adaptive
  selection/reordering uses explicit native maps. Ordinary fitting is unchanged.
- New-measurement projection estimates signed feature offsets for affine
  sources using fixed activities and an augmented ones column. Intercept-only
  sources are supported; dependent columns and excessive allocations are refused.
- Recognize audited generalized-binary and nonnegative-unimodal EBNM support;
  zero-factor native interval extraction returns a typed empty table.
- Add `fit_seminonnegative_ebmf()` with semantic sample/feature priors, native
  nonnegative activity and signed effects, both input orientations, and explicit
  support/refit handling for custom effect priors.
- Projected fits now use their own native factor indices, scaling, PVE,
  posterior SDs and sampler coordinates. Intersection projection reports
  unmeasured-feature uncertainty as unavailable, and returned native fits retain
  external offsets through standardization and serialization.
- Metadata tuples remain distinct when labels contain punctuation; the same
  identity rule protects donor bootstrap strata.
- Automatic stability fitting now requires an explicit callback for edited,
  projected or custom fits. Ordinary fits and pure reorders retain their recipes.
- Validate native factor maps and return typed empty summaries, rankings,
  specificity and marker correlations for valid zero-factor fits.
- `reorder_factors()` now returns a fit with a new identity, so a view made before
  the reorder is refused instead of silently pairing posterior SDs with the wrong
  factor. Provenance `engine` and `engine_version` are single strings (a partial
  `$` match had turned them into lists).
- Fix `fix_factors(ids = )`: partial IDs now fix exactly the named samples or
  features in both orientations (a single ID previously fixed the whole axis and
  several IDs errored). `fix_factors()` and `remove_factors()` now require
  `factors` instead of treating `NULL` as every factor.
- Tutorials state where their data and objects come from: a visible "Data and objects"
  section with the command to obtain the pancreas data, explicit reads of the cached
  files, a visible simulation in the protein article, and a README table of data
  sources. `restore-programs.R` now takes the cache as an argument.
- Every tutorial ends with a "Function credits" table built from `provenance()`,
  listing the functions it calls and their original authors.
- Internal simplification, no user-visible change: one engine check for optional
  packages, `.data` imported from ggplot2 instead of rlang, and the pancreas
  restore code shared across tutorials.
- Add tutorial coverage for every exported function: a new "Reorder, remove, fix
  and refit factors" article (`reorder_factors`, `remove_factors`, `fix_factors`,
  `refit_factors`, `factor_pve`) and examples for `fit_nonnegative_ebmf` (GTEx),
  `rank_features`, `factor_distinctiveness`, `rank_factors_by_metadata` and
  `factor_cor` (pancreas programs).
- `factor_gate_alignment()` gains `direction = c("high", "low")`; `"low"` compares
  a gate with the lowest-activity samples. The default is unchanged and the
  summary gains a `direction` column.
- Document that the `eligible_above = 0.5` default of `estimate_marker_thresholds()`
  is assay-specific (function manual and protein tutorial); no behaviour change.


## Backend compatibility (FU-0, 2026-10-18)

* `fit_ebmf()` now preflights the installed flashier's native capabilities at
  fit entry and fails fast with an actionable message (naming the pinned
  upstream revision) when a symbol the pipeline needs late is absent —
  previously a flashier lacking `flash_fit_get_p2m`/`flash_factors_remove`
  failed only after the expensive fit + nullcheck completed. Compatibility
  audit at the pinned runtime (flashier 1.0.60 @ fd44811): all 19 flashier
  symbols used by the package are exported. Trivial extractors are
  deliberately not preflighted.

## Grouped specificity validation (FU-2B, 2026-10-18)

* `factor_specificity()` validated without API change (return columns
  pinned by test): NA group/unit labels fail explicitly through the
  metadata guard ("Missing metadata values for IDs: ..."), never an
  `any(NA)` exception; AUC matches a hand-computed Mann-Whitney oracle and
  pROC in both directions; tied optimal thresholds are counted
  (`threshold_ties`) with the smallest finite threshold chosen; units
  spanning both classes are visibly double-counted
  (`n_positive_units + n_negative_units > n_units`) so non-independence is
  transparent; extreme imbalance (1 vs 19) stays typed and hand-checkable.
  Documentation now states explicitly that AUC/threshold/direction are
  in-sample descriptive statistics (never predictive, no CIs by default)
  and that unit counts are descriptive tallies when donors span groups.

## Membership burden, centering and global feature selection (FU-2A/FU-2C, 2026-10-19)

* `factor_diagnostics()` gains opt-in `membership_burden = TRUE`: per-sample
  b_i = sum_k max(A_ik, 0) / m_k over the raw resolved activities, with
  m_k = max_i A_ik; zero-max factors contribute 0 and are recorded
  (`zero_max_factors`), and the definition, maxima, clip policy and reference
  population are declared in the analysis metadata. Off by default.
* New `center_factor_activity(fit, metadata, by, baseline)` centers grouped
  activity means against an explicitly declared baseline
  (`equal_group` or `observation_weighted`) for display; consumes
  `summarize_factor_activity()`, records the baseline in the metadata, and
  attaches no inferential label.
* New `factor_features_global(fit, factors, n)` selects features globally
  across selected factors by max_k |effect| with stable ties (score, then
  feature ID), returning the argmax factor and signed estimate;
  `plot_factor_features(select = "global")` composes it.

## RNA-to-ADT fixed-activity projection validation (FU-3A, 2026-10-19)

* `project_factor_features()` validated for the CITE-seq configuration
  (Y_ADT ~ A_RNA B_ADT', activities fixed, exact cell intersection, signed
  protein priors, declared ADT units) against a native oracle: per-protein
  effect correlation > 0.98, held-out-cell reconstruction from RNA
  activities alone, explicit refusal of rank-deficient activities, of query
  cells absent from the reference, and of undeclared units. Worked example
  added to the protein-interpretation vignette. No numeric parity with
  ImmGenT-GP is claimed (different backfit controls).

## Chunked projection decision (FU-3B, 2026-10-19)

* Feasibility audit (reports/FU3B_CP6_CHUNKING_DECISION.md): exact
  row-chunked projection is a NO-GO under the current API because activity
  priors are re-estimated by empirical Bayes over the query rows (and
  `var_type = 0/2` pools variance over rows); exact chunking would require
  frozen fitted priors, an upstream flashier feature. Approximate chunking
  is not implemented; oversize queries keep failing fast against
  `dense_budget_bytes`. `project_factor_activity()` documentation now states
  the query-set dependence explicitly.


# flashier.utils 1.1.0

- Add signed program networks, annotation comparisons, descriptive diagnostics,
  native ROC specificity, protein thresholds and equal-population gates.
- Add native fixed-activity and fixed-program EBMF projection with explicit
  units, preprocessing, coverage, rank and allocation checks.
- Add workflow tutorials and a published-data pancreas example adapted from
  Peter Carbonetto's single-cell Jamboree analysis. Tutorials now use real data
  (flashier's GTEx matrix; pancreas caches from `.github/scripts/prepare-pancreas.R`,
  built by the pkgdown workflow and gated off during package checks).

# flashier.utils 1.0.0

* Added signed decoupleR source/target/mor networks through pathways, with native
  ULM/MLM scoring, no default centering, explicit one-to-one feature mapping,
  network hashes and retained annotation losses. Duplicate edges and invalid
  weights are refused; weighted networks cannot use fgsea. Scores annotate
  estimated feature programs and condition on the EBMF fit.
* Implemented the 34 planned exports around public flashier and optional engines.
* Enforced named coordinates, biological units and fixed-basis backprojection.
* Added donor resampling/checkpoints, independent-fit matching and leakage-safe
  projection/masked scoring with retained failures and explicit denominators.
* Added Seurat/uwot and display-only ECDF; original inferential scores stay signed.
* Added source receipts, complete helper notices and actual runtime records.
* Explicitly narrowed unsupported NNLM/cycling, SuSiE, CorShrink dosym/bootstrap,
  frozen-noise projection and brms/uncertain-response modes. No fallback engines.
* Confirmed package identity flashier.utils and maintainer Mikhael Manurung.
  GitHub distribution with Linux release/devel, macOS and Windows checks.
* Full-engine CI verifies every suggested dependency loads; installs GSL for mashr.
