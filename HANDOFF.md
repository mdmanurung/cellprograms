# Handoff — cellprograms: cross-cell-type sharing, adjustment, review fixes

**Date:** 2026-10-08 · **Branch:** fix/review-defects @ 2c7f980 (off main 5300e92; unpushed) · **Status:** PR #1 merged; review-defect fixes committed locally; benchmark work sits on unpushed `feat/sofa-comparator` (5 commits, tip 7efc410)

## Goal
Package-level, defensible cross-cell-type claim from per-cell-type EBMF programs (COMBAT first). Done = stratified/blocked sharing null, a multi-cell-type program object, and the cheap review fixes merged, with COMBAT conclusions rerun under them.

## Next action
`git push -u origin fix/review-defects` and open a PR (user's call), then add a `strata=` argument to `principal_angles`/`sharing_spectrum` (`R/principal-angles.R`: permute donors of b within strata inside `.pa_engine`, line ~66).

## State
- Uncommitted: `?? benchmarks/biomni_replication/results/sim_adj4/` (belongs to `feat/sofa-comparator`); nothing else
- `feat/sofa-comparator`: fdfc38a SOFA harness, 7f86492 SOFA vs arms (20 reps), 3027839 / 0766a0d / 7efc410 before-after recovery, ordinal severity, COVID-only baseline. Not pushed, no PR.
- Tests: full suite passes (covariates 14, review-fixes 13, sharing-calibration 13, ebmf 16, data-model 19, extract-stability 7, mofacellular 13, principal-angles 4). Run: `R_LIBS_USER=/nonexistent /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin/Rscript -e 'library(stats);library(utils);for(f in list.files("R",full.names=TRUE)) source(f); library(testthat); source("tests/testthat/helper-toy.R"); test_dir("tests/testthat", env=globalenv(), reporter="summary")'` (devtools fails: rlang too old)
- Running: nothing of mine. Use SLURM `--partition=all`. SOFA env: `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/sofa-env` (py3.8; set `NUMBA_CACHE_DIR` per task).

## Results (COMBAT, K=10, scores space; `benchmarks/biomni_replication/results/sharing/`)
- Median excess_frac by adjustment: none 0.50, Institute 0.49, Institute+Outcome(7 lvl) 0.36, Institute+ordinal severity+Sepsis indicator 0.43; COVID-only none 0.49 vs Institute+ordinal 0.49. All callable pairs stay shared; every non-call is underpowered.
- The 0.50 -> 0.36 drop is mostly the 6-df cost of `Outcome` in small cell types (GDT/DP/DC/PB), not severity biology: with 1-df ordinal, program loss = baseline instability (12 vs 11 lost of 97).
- Simulation (20 reps): unadjusted test rejects 100% with no sharing; hard and soft adjustment reject 0.15/0.00 and keep power; soft == hard (expected: near-flat prior); SOFA guided leaves batch corr 0.16, AUROC 0.85-0.95, unguided ~0.5. Scores not like-for-like.

## Locked decisions
- Hard adjustment (`covariate_mode="residualize"`) is the default; soft kept as option only (user, 2026-10-08)
- First COMBAT adjustment technical only (Institute); then Institute+Outcome; then ordinal severity (user, 2026-10-08)
- Do not port L1 shared priors, L2 init borrowing, EV-BIDIFAC, ICA/PVA; do NOT replace repo `fit_celltype_programs` (user, `biomni/REVIEW.md`)
- Repo defaults K=10, `point_laplace`, `var_type=1` (user, 2026-10-07)
- Partition `all`, not `medium` (QOS caps 4 cpu/8G)

## Dead ends — do not redo
- `p_pair` alone to compare COMBAT runs: saturates at the permutation floor for 41-44 of 45 pairs; use `excess_frac`/`z_pair`.
- `ebnm_flat` for fixed covariate factors: NA ELBO in flashier; use fixed-scale `ebnm_normal`.
- AUROC of z_pair as evidence: 1.00 even unadjusted where FPR = 1.0 (Biomni D7); report FPR/power.
- 7-level `Outcome` as primary adjustment: fuses disease type and severity (blank = 23 Sepsis donors), 6 df in ~50-donor cell types.
- Bootstrap with replacement for stability; fixed seed inside a data-generating loop (resets RNG).
- Raw-count pseudobulk, z-scored baselines, Biomni 10-dim Spearman floor (see `biomni/REVIEW_round2.md`).

## Read first
- `biomni/REVIEW_round2.md`; Biomni reviews of the sharing test (tsk_015R57NqiKS1IVa2f00lOlJX) and of the implementation (tsk_0165grpXiyKQ1bHR90MnLIvG, defects D1-D20, in the Biomni web app only)
- `R/principal-angles.R` (`.pa_engine`, `sharing_spectrum`), `R/ebmf.R` (`.cov_design`, fallback ladder at ~line 143)

## Open items
- Not fixed from the review: stratified/blocked null (D1, next action); multi-cell-type program object / program graph (D12); cell-count weighting via flashier `S` (D5); common gene panel and post-adjustment gene selection (D6); K cap <= n-q-1 (D7c); optimal (not greedy) matching (D13a); zero-fill downstream (D4, available-case).
- K=10 for every cell type makes small types uncallable; per-cell-type K by split-half stability. DC COMBAT fit used `point_normal+no_nullcheck`; treat DC pairs as sensitivity.
- Earlier COMBAT p-values used one seed for all pairs; rerun with per-pair seeds before quoting.
- Biomni project also has an unrelated task "Review EGSV2 P vs NP findings F-003"; ignore.
