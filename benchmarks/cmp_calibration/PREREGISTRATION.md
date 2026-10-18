# CMP calibration audit — preregistration (CP-3A, spec v2 family A-D)
# Date: 2026-10-18. Frozen BEFORE any confirmation run. Development smoke
# (20 replicates/scenario, n_perm = 199) may only catch code failures and
# size the confirmation arm; no calibration claim is made from smoke runs.

## Design

Score-level simulation driving the PRODUCTION find_coordinated_programs()
entry point on synthetic `fit$scores` (the EBMF layer is upstream of the
statistical question: calibration lives in rank-normalization, block
residualization, pairwise-complete correlation, the node max-statistic and
the within-block permutation null).

- alpha = 0.05, tau = 2, link_threshold = 3, min_shared = 20, min_donors = 20
  (package defaults; min_shared/min_donors lowered to 10 for the sentinel
  design with N = 30 donors).
- Sentinel design: C = 2 cell types, K = 4 programs/type, N = 30 donors.
- Realistic design: C = 6 cell types, K = 8 programs/type, N = 100 donors.
- Confirmation n_perm = 499 (p-value resolution 1/500 = 0.002 << alpha).

## Scenarios

- A0 global null, iid: within-type equicorrelation 0.3, no cross-type
  dependence, single block. EXPECTATION: FWER ~= alpha.
- A1 global null, blocked: 3 equal sites, site mean shifts only (valid
  exchangeability within blocks; residualization removes block means).
  EXPECTATION: FWER ~= alpha.
- B1 partial null: one planted cross-type driver (rho = 0.6) linking one
  program in each of 3 cell types; all other programs independent.
  METRICS: null-node FPR (strong-FWER challenge) and planted-node power.
- C1 block equicorrelation (exchangeability violation, D-1 replication):
  donor-level block component shared ACROSS cell types within blocks
  (rho_b = 0.55), unremoved by mean residualization. EXPECTATION:
  anti-conservative; quantified, not required to pass.
- C2 shared donor-dependent missingness: top-quartile latent donors missing
  in ALL cell types simultaneously (informative, shared across types).
- C3 unequal blocks + singleton block (sites of size 40/8/1/1).
- C4 zero-variance program injected (constant scores) under global null.
- D1 grouping: planted 3-type CMP with a signed member (driver signs
  +1/+1/-1). METRICS: membership precision/recall/Jaccard, sign recovery vs
  planted.
- D2 grouping stress: three independent pairwise drivers forming a
  sign-inconsistent triangle (A-B +, B-C +, A-C -; sign product -1).
  METRICS: whether grouping/orientation reports a coherent consensus for an
  incoherent graph (it must not market an arbitrary sign as inferred
  consensus); membership ambiguity recorded descriptively.

## Decision rules (frozen)

- Report empirical FWER/FPR with exact binomial 95% Monte Carlo CIs.
- A0/A1 support the global-null calibration claim if the CI covers alpha.
- No strong-FWER claim is made regardless of outcome unless B1 null-node FPR
  CI lies at/below alpha; otherwise documentation narrows the claim (this is
  a documentation outcome, not a code-change trigger).
- C-family results are reported as characterized violations with quantified
  anti-conservatism; they gate WARNINGS/documentation, not silent patches.
- D1 reports descriptive grouping accuracy; no inferential claim on groups.
- Confirmation replicate count: 500 per scenario (SE ~ 0.0098 at p = 0.05;
  increased to 1000 if measured runtime permits within the compute budget).
- Confirmation seeds: 202610180 + scenario index * 1000 + replicate.
  Development smoke seeds: 9000 + replicate (disjoint).

## Budget

Smoke: 20 replicates x 8 scenario-design combinations, n_perm = 199.
Confirmation: sized from measured smoke runtime; target <= 4 h wall on the
session worker, parallelized over replicate batches. If the budget cannot
hold 500 replicates/scenario, the confirmation arm is downscoped and the
shortfall is reported explicitly rather than claimed as passed.
