# 04 — Experiment Plan: M3 Battery, Approach-B Experiments, M4-M6 Scaling

Workstream G / pre-registered BEFORE any M3 run executes (commit timestamp is the
registration timestamp). Elaborates the approved brief plan; companion to
`02-statistical-design.md` (criteria A1-A5, B1-B3, C1-C4; leakage rules B-L1..L3) and
`03-computational-design.md` (cost model, equivalence map, checkpointing policy). No
registered criterion is relaxed; the M3 gate remains "demonstrated recovery at reduced
cost, or documented failure".

Delegated decisions recorded as assumptions (user delegation, 2026-10-17): (1) scope =
full battery; (2) this document is committed before any M3 run; (3) seed budget = 10
seeds synthetic / 5 seeds real, >=3 repeats for timing claims; (4) cNMF is added as a
reference anchor (S5) in the recovery benchmark.

## 1. E-M3.1 — Initialization-strategy benchmark (realizes approach C)

**Synthetic generator** (new Python generator; the workstream-E generator, 80x600 and
donor-less, is too small for cell-level strategies — same statistical design, scaled and
donor-structured):
- 5,000 cells x 1,500 genes; 50 donors x 100 cells; donor ids recorded (donor = unit of
  dependence, B-L2).
- Truth: K_true = 3 broad programs (150 marker genes each, block-structured W; donor-level
  usage z with donor random effect sd 0.3) + 1 rare program (5% of cells, 80 markers) +
  2 null factors (pure noise columns).
- Counts = round(exp(Y)), Y Gaussian latent as in workstream E; nested Bernoulli dropout
  at rates 0 / 0.1 / 0.3 (u drawn once per rep; nested masks).
- Donor covariate z_age (standardized, donor-level); planted moderation on program-1
  usage amplitude, beta in {0, 0.3, 0.6, 1.0}; E-M3.1 runs on beta = 0 data, E-M3.2 on
  the beta grid. Fixed seed stream recorded in the result contract.
- Synthetic pseudobulk = donor-level aggregation (no cell types planted); metacells =
  k-means within donor on normalized expression, n_cells weights (workstream-F convention).

**Arms** (matched preprocessing = registered shifted-log recipe; operating point K=10,
laplace/laplace, COLUMN_WISE, backfitting + pruning, seed family 20261016):
- S1 from-scratch single-cell fit (reference).
- S2 pseudobulk-init + chunked project to all cells.
- S3 metacell-init + chunked project.
- S4 init + refine: user-init passthrough from S2/S3 factors, full EM continue (maxit 50,
  ELBO-plateau stop).
- S5 cNMF anchor: k search 2-8 + consensus, same matrices, zero-gene drop recorded.

**Donor-safe evaluation protocol (binding, all arms)**: fit on the 80% training donors
only; held-out prediction = project the 20% test donors' cells with F fixed, evaluate
Pearson cor + RMSE on a fixed random 10% entry mask of the test donors' shifted-log
matrix; recovery metrics computed on training cells (synthetic truth) or vs S1 (real).

**Metrics per arm** (recorded in the result contract): Hungarian-matched |cor| of scores;
top-50 gene Jaccard per matched program; held-out prediction (above); seed stability
(pairwise matched |cor| across seeds); wall-clock; peak RSS; ZeroRate; K_eff after
pruning (C4: no discovery claim without ZeroRate + K_eff).

**Gates (registered in 02 §C)**: C1 — S4 vs S1: ELBO within 1% at >=5x lower wall-clock
at matched maxit; C2 — rare program |cor| >= 0.7 in >=8/10 seeds via C1->C2->C3 (S1
reference rate reported either way); C3 — cross-granularity concordance >= the flashier
workstream-F pilot (same metrics, same data); C4 — transparency as above.

**Real-data arms**: UCSF CD8 granularities from workstream F (pseudobulk 99x2000;
metacell 2138x2000; singlecell 29,998x1,934 after zero-gene drop), 5 seeds/arm. The
real-data cNMF anchor is REUSED from workstream F (no duplicate runs). S1 feasibility
measured in M2: 237 s at 10k cells, ~linear in N (~12 min at 24k train cells, ~2.6 GB).

## 2. E-M3.2 — Approach-B powered covariate-moderated experiments

- Grid: beta in {0.3, 0.6, 1.0} x n_donors in {40, 99} x 10 seeds (beta = 0 arms are
  controls; generator as §1).
- Arms: B-fit (lcash prior on L, X_l = z_age), A-fit (covariate-free), A + post-hoc
  regression of scores on z_age (the B3 comparator).
- Evaluation on held-out donors only (B-L3). B1: B-fit held-out Spearman(usage, z_age)
  >= 0.8 on the moderated program while A-fit |Spearman| <= 0.5. B2: non-moderated
  programs still meet A1 thresholds under the B-fit. B3: if post-hoc regression matches
  moderation, the report says so.
- Outputs: smallest detectable beta (power curve); whether prior moderation beats
  post-hoc regression (the paper's value proposition, tested not assumed). The completed
  age pilot (negative at the beta~0 operating point) is cited as motivating evidence.
- Real-data arm: donor age, 99 UCSF donors, 5x2-fold donor splits, labeled pilot-grade.

## 3. E-M3.3 — Scaling (M4-M6) and E-M3.4 — R contract test

- **M4**: 1M x 3000 full chunked path (discovery on 50k subsample + projection of 1M),
  background, stage checkpoints per 03 §6 (integrity-verified before claimed). Estimate
  from M2 measurements: projection ~26 min (157 s/100k, 3 sweeps), discovery ~20 min,
  total < 2 h, memory flat ~4-6 GB.
- **M5**: 500k x 3000 shard profile (5x the M2 measurement) to confirm linear time / flat
  memory; final 5M extrapolation table includes the 5M-dense = infeasible-on-64GB row.
- **M6 gate review**: adopt/revise/abandon recommendation with all evidence tables
  (E-M3.1, E-M3.2, M4, M5), upstream contribution drafts (silent-divergence report, 8/30
  grid; frozen-prior feature request), and the Phase F benchmark handoff.
- **E-M3.4** (plan test 6): Python fit -> R cell_program_fit -> canonicalize_programs on
  worker-1 (R 4.4.2); runs once during M3.

## 4. Compute, execution targets, deliverables, acceptance

- Estimates (measured bases): E-M3.1 synthetic ~6-8 h single-worker -> ~3-4 h across
  worker-0 + worker-3 (adapter staged via shared mount, PYTHONPATH=src); real arms ~1-2 h;
  E-M3.2 ~2-3 h (~180 fits x 30-60 s); M4 ~2 h; M5 ~0.5 h; E-M3.4 ~0.5 h. Total ~9-13 h,
  CPU-only, background with checkpoints; no HPC.
- Execution targets: worker-0 (dev tree, torch 2.7.1+cpu), worker-3 (staged adapter copy),
  worker-1 (R contract test; canonical repo for commits).
- Deliverables: this document (first commit); per-arm CSVs + PNG figures to the results
  store and shared workspace (md5-verified both sides); ledger entry per experiment;
  commits per milestone; final = M6 gate review with the adopt/revise/abandon
  recommendation.
- Acceptance: fixed seeds, pinned revisions, recorded backend versions; every checkpoint
  verified by enumeration + size before reported; criteria evaluated exactly as
  registered — pass/fail recorded, no post-hoc relaxation; failures documented as
  documented-limitation findings feeding the recommendation.
- Workstreams E (aggregation pending) and F (two arms running) continue autonomously and
  fold in as pilot evidence without blocking this plan.
