# cellprograms — limitation register (REL-1, 2026-10-19)

Scope: user-visible limitations established by the spec-v2 milestone evidence.
Each entry names the evidence artifact. This register is normative for
documentation and claims; do not remove entries without a confirming
re-analysis.

## L1. Projection is query-set dependent (CP-6 / FU-3B)

`project_programs_to_cells()` delegates to
`flashier.utils::project_factor_activity()`, which re-estimates activity
priors by empirical Bayes over the query rows (and pools residual variance
over rows at `var_type = 0/2`). The same cell projected within two
different query sets can receive different activities. Chunked projection
is an approximation, not an exact decomposition, and is not offered;
oversize queries are refused pre-flight against `dense_budget_bytes`.
Evidence: flashier-utils `reports/FU3B_CP6_CHUNKING_DECISION.md`.

## L2. CMP context association: C-family FWER at the boundary (CP-3A/CP-3B)

Confirm-scale simulation (5500 replicates, n_perm = 499): A0/A1 nulls are
calibrated (FWER 0.044–0.054, CIs cover 0.05); B1 supports strong FWER
control (null-node FPR 0.001 [0.001, 0.002]); the C-family (correlated
predictors) sits at FWER 0.040–0.058, CIs covering or near alpha — a
boundary result that is documented rather than patched, per preregistration
(CP-3B is docs-only). D1 plant power 1.000; D2 0.825 [0.805, 0.844].
Evidence: `benchmarks/cmp_calibration/CONFIRM_ANALYSIS.md`.

## L3. Program sharing: small-support cross-group programs are not promoted (CP-4)

At confirm scale (20 reps/arm), planted sex/genotype-shared programs
(N3/N4) were never promoted as shared in any arm — n_top/support dilution
(Amendment 1). The sharing detector is conservative for small-support
cross-group programs; absence of promotion is not evidence of absence.
Posthoc artifact removal partially promotes ambient contamination (0.600
[0.361, 0.809]) — descriptive only. All preregistered gates PASS.
Evidence: `benchmarks/program_sharing_sensitivity/CONFIRM_ANALYSIS.md`.

## L4. Frozen preprocessing specs do not carry measurement units (CP-5)

`freeze_preprocessing()` records feature sets, design and centering — not
`feature_scale`. A refit under `preprocessing_spec` therefore keeps
`feature_scale = "unknown"`, and downstream projection from such a fit is
refused (units undeclared). Declare units at discovery-fit time if the fit
will serve as a projection reference. Evidence: CP-5 vignette build.

## L5. Cross-cohort sign/scale comparability (CP-5)

Program sign and scale are comparable only within a transfer pathway
(fixed-basis projection, or de-novo fit after matching), never across
pathways or across de-novo fits. The replication rule must be
same-sign-as-discovery, not a fixed direction.
Evidence: `vignettes/articles/cross_cohort_replication.Rmd`.

## L6. Metacells: first-use environment build (fold-in)

`MetacellsByGroupsSEACells()` uses a basilisk-managed Python environment
registered under pkgname `cellprograms`; the first call after install
rebuilds the environment (one-time cost). mcRigor is a GitHub Suggests
dependency (JSB-UCLA/mcRigor), not CRAN. Evidence: fold-in ledger entries
2026-10-18.
