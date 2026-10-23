# flashier.utils — limitation register (REL-1, 2026-10-19)

Scope: user-visible limitations established by the spec-v2 milestone
evidence. Normative for documentation and claims.

## L1. Projection is query-set dependent; no exact chunking (FU-3B)

`project_factor_activity()` / `project_factor_features()` re-estimate
priors by empirical Bayes over the query set; `var_type = 0/2` additionally
pools residual variance over rows. The same query row in two different
query sets can receive different estimates. Exact row-chunking would
require frozen fitted priors (upstream flashier feature); approximate
chunking is not implemented. Oversize queries fail fast against
`dense_budget_bytes`. Evidence: `reports/FU3B_CP6_CHUNKING_DECISION.md`.

## L2. RNA-to-ADT parity is a native-oracle validation (FU-3A)

The CITE-seq configuration (fixed RNA activities, signed protein priors,
declared ADT units) is validated against the native fixed-activity
projection itself (per-protein effect correlation > 0.98; held-out-cell
reconstruction from RNA activities alone). No exact numeric parity with
ImmGenT-GP is claimed: backfit controls and iteration schedules differ.
Evidence: `tests/testthat/test-fu3a-rna-adt.R`.

## L3. Specificity and burden statistics are in-sample descriptive (FU-2A/FU-2B)

`factor_specificity()` AUC/threshold/direction are in-sample descriptive
statistics — never predictive, no confidence intervals by default; donor
units spanning both classes are visibly double-counted
(`n_positive_units + n_negative_units > n_units`). `factor_diagnostics(
membership_burden = TRUE)` is a descriptive raw-activity diagnostic;
zero-max factors contribute 0 and are recorded. `center_factor_activity()`
is display centering against a declared baseline with no inferential
content. Evidence: `tests/testthat/test-fu2a-fu2c.R`,
`tests/testthat/test-fu2b-specificity.R`.

## L4. Unsupported modes remain narrowed (FU-0 and earlier)

NNLM/cycling, SuSiE, CorShrink dosym/bootstrap, frozen-noise projection and
brms/uncertain-response modes are explicitly refused; there are no fallback
engines. Backend compatibility with the pinned flashier (1.0.60 @ fd44811)
is preflighted at fit entry. Evidence: NEWS 1.1.0.9000, FU-0 audit.
