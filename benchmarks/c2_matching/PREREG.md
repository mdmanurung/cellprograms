# PREREG: C2 registered check — optimal vs greedy program matching on stability refits

Status: **DRAFT** (registration is a scientific commitment; submitted for user
review before the main run — plan assumption: rule-determined outcomes
auto-proceed, this draft documents the rule before the check runs).

## Question

`match_programs()` matches programs between fits greedily
(`.greedy_match()`: repeatedly take the global-max |cosine| pair, remove the
row and column). Greedy assignment can be suboptimal for the total matched
weight. Does optimal (maximum-weight bipartite) matching change the matched
set in practice?

## Registered decision rule (ROADMAP M6)

Adopt optimal matching into `R/` (replace the greedy step in
`match_programs()`) **iff** the pooled discordance rate is **>= 5%** of
programs. Otherwise record "dropped; rule did not fire" in HANDOFF/ROADMAP
and keep greedy matching.

## Data

The Stage-A real CMP-detection fits (cached, not committed):

- `stephenson_pb_sex_K30_fit.rds` (12 cell types, ~100 donors/ct, 2000 genes,
  K = 27-30, covariate_mode = residualize)
- `combat_pb_sex_K30_fit.rds`

## Refit protocol

Exactly `assess_program_stability()`'s procedure, per fit:

- n_boot = 10 subsample refits per fit (pre-registered extension: +10 boots
  iff the pooled rate lands in [0.02, 0.10), to separate a near-threshold
  rate from noise; the extension is decided on the pooled rate of the first
  10 boots before any unblinding of per-cell-type detail).
- Per boot: subsample donors per cell type at `resample_frac = 0.8`
  (minimum 5 donors), refit with the fit's stored settings
  (`seed = provenance$seed + b`), canonicalize.
- Deviation from `assess_program_stability()`'s single RNG stream: each boot
  seeds its own donor subsample with `set.seed(b)` so boots are independent
  and the run is resumable per boot. The fit seed convention is unchanged.
- Refits that error are dropped (as in `assess_program_stability`); the
  denominator counts only completed (boot, cell type) pairs.

## Matching comparison

Per (boot, cell type) with K > 0 on the original side:

- Similarity: |cosine| between original and refit loadings on the intersected
  gene panel — the same matrix `match_programs()` computes.
- Greedy assignment: `.greedy_match()` (the package's current matcher).
- Optimal assignment: maximum-weight bipartite matching on the same matrix
  (`clue::solve_LSAP`, transposed when the refit side is larger), i.e. the
  assignment maximizing total |cosine|.
- Threshold: keep pairs with |cosine| >= 0.5 (package default
  `sim_threshold`), applied after assignment in both methods.
- A program is **discordant** if its post-threshold partner on the refit side
  differs between greedy and optimal (including one method leaving it
  unmatched).

## Metric

Pooled discordance rate = discordant programs / all original-side programs,
over all completed (boot, cell type) pairs of both fits. Secondary: total
matched weight under each method (how suboptimal greedy is in total |cosine|).

## Cost

Pilot-measured; see run_check.R header and HANDOFF for the measured
per-refit time and the sandbox vs Slurm decision.
