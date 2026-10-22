# prereg-coord-v1 — Cross-cell-type coordination benchmark (C1)

**Tag:** `prereg-coord-v1` · **Date:** 2026-10-21 · **Status:** FILLED — all margins frozen from the runtime and amplitude pilots before the main run
**Supersedes:** the informal 3-arm nomination (2026-10-21) and the v2 candidate plan
(`coordination_benchmark_plan_v2.md`), whose mandatory changes (P0-1 endpoint redesign, P0-3
priv_corr mode, scorer harness, adjustment position, FPR gates) are incorporated here.

## 1. Scientific question

Given per-cell-type donor-level pseudobulk expression (genes × donors, per cell type) with
donor-level nuisance structure, which integration strategy recovers **coordinated programs** —
latent programs whose donor activity is shared across two or more cell types — with the best
membership recall/precision at controlled spurious-call rates?

The decision this informs: whether the CellPrograms package should integrate per-cell-type EBMF
fits post hoc (pairwise scores-space tests), jointly at expression level (block EBMF), or jointly
with an empirical-Bayes tensor factorization (EB-CP), before real-cohort coordination claims are
made.

## 2. Arms

| # | Arm | Class | Role |
|---|-----|-------|------|
| 1 | `pa_cc` | pairwise, scores space, principal-angles sharing spectrum + connected components | baseline (decision-rule anchor) |
| 2 | `block_E` | joint, expression level, block EBMF (shared-first ordering), block-mass membership | primary challenger |
| 3 | `eb_cp` | joint, expression level, EB tensor CP (greedy rank-1 + backfit), model-based sparse membership | gated challenger — **gate passed** (§3) |
| A1 | `subspace_sum` | joint, scores space, multi-way eigenproblem | attribution: joint-vs-pairwise at scores space |
| A2 | `gene_match` | pairwise, loadings space, cosine + gene-label permutation | attribution: gene-usage mechanism |
| A3 | `node_cor` | pairwise, scores space, correlation + max-T FWER | attribution: simple calibrated correlation (post Defect-D fix) |
| — | `oracle` | true-loadings projection + node_cor machinery | detectability floor, **node-based arms only** (fix 8: no joint-arm oracle is defined) |

**Excluded, with reasons:** CMP on main (four P0 audit findings, CP-3A not started, no calibration
evidence); PMD (uncalibrated permutation null, rejects everything); evb_E (convergence sensitivity —
max_iter ≥ 1000, 0.5–75 min, module splitting at default budget — and no real-data expression-level
evidence); `pa_prog`/`score_cor` (never implemented; superseded by this arm list).

## 3. eb_cp gate record (P0-2)

Gate criteria from the v2 plan and their disposition:

- (a) Implemented: flash-style greedy rank-1 CP with matricized EBNM normal-means updates per mode
  (point-laplace on gene and cell-type modes, normal on donor mode), SVD-seeded, warm-started
  full-data refit + 2 backfit sweeps. No monotone joint-ELBO claim is made; selection is by
  held-out likelihood, not training ELBO.
- (b) Smoke-tested on `gen_coord` null + act+genes families (`test_eb_cp.R` §1, §3–4).
- (c) Null behavior: **1/20 null datasets produced any coordination call** (gate: ≤ 1/20).
  Borderline-pass; the main run measures the actual spurious rate on 16 null replicates (§6).
- (d) Membership rule frozen: component member = cell type with |w_ct| ≥ 0.2 · max|w| over the
  cell-type mode; a coordination call requires ≥ 2 members.

Component-selection rule (frozen): donor-synchronous 50/50 split, **stratified on institute**;
a greedy component is kept iff its held-out association z > 3, with B-side per-row
self-studentization (each held-out row's own marginal variance). The stratified split is
load-bearing: an unstratified split leaves residual institute imbalance in each half after
full-data residualization, and the greedy pass then consumes components on ct-specific nuisance
directions (probe: identical inputs, stratified → clean global R100 component; unstratified →
nuisance-fragmented trajectory, 0 calls). Variance model: per-(cell-type, donor) — the tensor
analogue of flashier `var_type = 1`, the package's validated pseudobulk default.

## 4. Data-generating process

Simulator: `gen_coord` (benchmarks/coordination/sim_coord.R). Frozen configuration:

- 10 cell types (ct01–ct10); 100 donors; per-ct donor counts 100 ×4, 75 ×3, 50 ×3, with a shared
  30-donor core among the 50-donor cell types.
- G = 2000 genes (frozen from the runtime pilot; also the simulator headroom floor — worst-case
  gene-pool draw is 1850 for the co-planted family: 850 program + up to 1000 private genes);
  50 program genes per planted program.
- Nuisance: institute (3 levels, probabilities .60/.25/.15), effect `1.5 · outer(bz, bload)`,
  sparse random gene loadings; per-donor heteroskedastic noise `sd_i = sqrt(median(nc_i)/nc_i)`.
- Private programs in every cell type (priv_amp = 1): truth treats them as distinct programs.
- `priv_corr` mode: 3 decoy pairs of private programs in distinct cell-type pairs with correlated
  donor activity (cor_priv = 0.7) and disjoint genes — the S5 trap. A multi-cell-type call over a
  decoy pair is a **trap failure**.
- Breadth families: `R100` (10 member cts), `R50` (5), `R20` (2), `R20small` (2, both small cts);
  `all` co-plants R100+R50+R20 (secondary/sanity only — pa_cc's connected-component output is
  degenerate under co-planting, P0-1).
- Signal modes: `act_only` (shared activity, disjoint genes — primary), `act+genes` (shared
  activity + shared genes), `genes_only` (shared genes, independent activity), `act+half`
  (secondary; first mode cut if over budget), `priv_corr`, `null`.
- Activity correlation across members: act_cor ∈ {1, 0.9, 0.7} sweep on the primary family
  (fix 7: act_cor < 1 violates the joint arms' rank-1 activity assumption; the sweep exposes the
  breakdown point and separates modeling failures from detection failures).

**Amplitude (per breadth, act_only):** chosen by oracle pilot — the smallest amp with oracle
pairwise power ≥ 0.5 under the frozen node_cor machinery (stratified permutation, max-T FWER,
α = 0.05), subject to stage-1 recovery ≥ 0.5 (fraction of (program, member-ct) cells with a
canonical stage-1 program whose scores correlate |r| ≥ 0.5 with true activity; on violation,
move one grid step up and record the deviation). Other modes reuse the act_only amplitude for
the same breadth.

**Pilot results (384 + 192 oracle datasets; 304 stage-1 fits, G = 2000, act_only):** oracle
pairwise power undergoes a sharp transition (R100: 0.06 → 0.40 → 0.80 at amp 0.10/0.15/0.20;
R50: 0.06/0.46/0.80; R20: 0.08/0.33/0.92; R20small: 0.00/0.33/0.46 → 0.96 at 0.30). Oracle
picks: 0.20 (R100/R50/R20), 0.30 (R20small). Stage-1 recovery is 0.00–0.12 at amp ≤ 0.30 and
0.85–0.88 at amp 0.50 — stage-1 estimation, not detection, is the binding constraint.
**Deviation recorded:** the stage-1 veto moves every breadth above its oracle pick.
**Final amplitude = 0.50 for all breadths** (oracle power 1.00; stage-1 recovery 0.85–0.88).
The race is thereby placed in the estimation-bound regime by design: the detectability floor is
at ceiling, and arms are discriminated on how much achievable signal their estimation and
integration pipeline preserves — the regime in which joint expression-level arms can exceed the
stage-1-imposed recall cap (~0.85) that bounds all stage-1-dependent arms.

## 5. Pipeline (frozen positions)

1. **Adjustment (fix 4):** residualize expression on institute **pre-stage-1**, identically for
   all arms. No arm receives fixed donor-mode covariate factors (the "both or neither" rule;
   neither). Residualizing fitted scores instead of data is excluded (stage-1 capacity is spent
   on the nuisance).
2. **Stage-1:** `fit_celltype_programs` — point-laplace loadings, max_factors = 10, var_type = 1,
   gene-centered, zero-variance genes dropped — then `canonicalize_programs`. Stage-1 `Y_used`
   matrices are bit-identical to standalone residualize+center output (verified 2026-10-21), so
   expression-level arms (block_E, eb_cp) and scores-space arms (all others) consume the same
   adjusted data.
3. **Arm-specific frozen settings:**
   - `pa_cc`: sharing_spectrum in scores space, stratified permutation, n_perm = 1999; calls =
     connected components over cell types of significant pairs (program_id = NA; activity = NULL —
     decoy evaluation not applicable, a preregistered limitation).
   - `node_cor`: per-pair score correlation, stratified permutation of row order on the a-side,
     max-T FWER over all pair × program edges, α = 0.05, n_perm = 999, min_shared = 20 donors.
   - `subspace_sum`: joint scores-space eigenproblem, J = 15, n_perm = 999, min_shared = 20.
   - `gene_match`: loading cosine over shared genes, gene-label permutation, max-T, n_perm = 999,
     min_genes = 100.
   - `block_E`: joint block EBMF at expression level, order = shared_first, max_factors = 20
     (non-binding: the runtime pilot's fit terminated at the numerically-zero factor 20),
     block-mass > 0.1 membership, call = block with ≥ 2 member cell types.
   - `eb_cp`: as §3 (Kmax = 15, heldout_z = 3, member_thresh = 0.2, stratified split, per-(ct,
     donor) variance, 2 backfit sweeps).
   - `oracle`: per-ct projection of adjusted matrices onto true loadings, then the node_cor
     machinery (n_perm = 999).

## 6. Endpoints and scoring

**Matching rule (frozen before pilot):** one-to-one Hungarian assignment (clue::solve_LSAP) between
emitted clusters and planted programs on membership-set F1 (padded square); a planted program is
**found** iff its matched cluster has recall ≥ 0.5 and precision ≥ 0.5; unmatched clusters are
**spurious**.

- **Primary endpoint:** membership recall, precision, and breadth error on **R50 and R20 in the
  each-alone family, act_only mode, pilot amplitude** (single primary cell family, fix 6).
- **Secondary:** R100 each-alone (sanity); act+genes and genes_only modes; R20small; co-planted
  `all` family (granularity caveat); act_cor sweep; runtime per arm.
- **FPR endpoints (fix 5):** spurious calls per emitted cluster **and** per dataset, where the
  per-dataset spurious rate = proportion of datasets with ≥ 1 spurious call (FWER-style; mean
  spurious calls/dataset also reported) on the null family (16 replicates); decoy rate (spurious
  cluster whose activity correlates |r| ≥ 0.5 with nuisance) and **trap rate** (call containing
  both members of a priv_corr decoy pair) on priv_corr (24 replicates).

## 7. Decision rule

Fixed hierarchy, evaluated in order: **eb_cp → block_E → attribution arms** (subspace_sum,
gene_match, node_cor) against the `pa_cc` baseline.

- **Adoption margin:** a challenger wins the primary endpoint iff paired membership recall exceeds
  pa_cc by **≥ +0.10** (project B1–B3 convention) at the pilot amplitude, with precision not worse
  by more than 0.05.
- **FPR gate (absolute + relative):** per-dataset spurious rate ≤ **0.10** AND ≤ pa_cc + 0.03.
  The absolute ceiling is binding even if pa_cc is degenerate (its documented estimated-nuisance
  behavior is FPR ≈ 1.0).
- **Oracle floor:** node-based arms' recall is interpreted against oracle pairwise power at the
  same amplitude (floor = what perfect stage-1 estimation would yield). No oracle is defined for
  the joint arms; their recall is interpreted absolutely.
- **act_cor interpretation:** a challenger that wins at act_cor = 1 but loses at 0.7 is a modeling
  loser (rank-1 activity assumption violated), not a detection loser; stated in advance so the
  sweep result cannot be rationalized either way.
- **Gene-sharing interpretation (same logic, second axis):** eb_cp fits one gene-mode vector per
  component shared across cell types; act_only plants disjoint gene sets per member, so the shared
  gene mode is diluted across non-member genes (worst for R100: each ct carries signal in 50 of
  500 union genes). A challenger that wins act+genes but loses act_only is a gene-sharing modeling
  loser. Runtime-pilot foreshadowing (one dataset, amp 0.5, act_only R50): eb_cp 0 calls vs
  pa_cc 1; and in a smoke replicate at amp 0.75 act_only R100, node_cor found the program
  (recall 0.7) while eb_cp emitted none — versus clean eb_cp recovery in act+genes at the same
  amplitude (gate test §3).

## 8. Replicate plan (Monte Carlo)

Shared dataset seeds across all arms within a replicate (paired comparisons).

| Cell family | Replicates |
|---|---|
| Primary: {R50, R20} × {act_only, act+genes} | 32 each |
| null | 16 |
| priv_corr | 24 |
| genes_only (R50) | 16 |
| R100 sanity (act_only) | 8 |
| R20small (act_only) | 16 |
| act_cor sweep {0.9, 0.7} × {R50, R20} act_only | 16 each |

Total = 272 datasets. Measured per-dataset cost at G = 2000 (runtime pilot, one act_only R50
dataset): stage-1 212 s, pa_cc 10 s (1999 perms), node_cor 4 s, subspace_sum 3 s, gene_match
30 s, block_E 593 s, eb_cp 9 s → ≈ 860 s/dataset; wall-clock ≈ 5.4 h on 12 parallel workers
(block_E is ≈ 70% of the budget). The primary family contributes 64 replicate–program pairs per
arm pair (32 reps × {R50, R20}); shared seeds make comparisons paired, removing dataset
variance; the worst-case Monte Carlo SE on the paired recall difference is 1/√64 ≈ 0.125
(realistic paired SDs are far below 1), and the confirmatory scale (96 reps → 192 pairs,
worst-case ≈ 0.072) tightens this below the 0.10 adoption margin. **Execution-scale note:**
this run is development-scale; the preregistered confirmatory scale for cluster execution is
96 replicates per primary cell (to be run if the development-scale decision is not already
decisive).

## 9. Known limitations (stated in advance)

- eb_cp null calibration is borderline (1/20); the main run's 16 null replicates are the real
  measurement, and eb_cp can fail the FPR gate even though it passed the milestone.
- pa_cc emits no program_id or activity: decoy/trap endpoints are not evaluable for it.
- Scores-space arms are known to degenerate at K ≈ 25–30 programs/ct with ≤ 50 shared donors
  (real-data audit); conclusions are scoped to the low-K regime except for one secondary high-K
  scenario if budget allows.
- Simulated institute nuisance is mean-shift only; real nuisance with interaction structure is
  out of scope for v1.
