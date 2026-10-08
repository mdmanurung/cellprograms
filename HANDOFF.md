# Handoff — Biomni review, replication, package follow-ups

**Date:** 2026-10-08 · **Branch:** feat/covariates-and-sharing-tools @ 44f20a6 (8 commits ahead of main b97fcad, not pushed) · **Status:** covariates, calibrated sharing test, effect size, stability/extract all in the package and tested; first adjusted-COMBAT result in; nothing from `benchmarks/` or `biomni/` committed

## Goal
Decide what from Biomni's `cellprograms` fork to bring into this repo, backed by a replication that can be trusted. Done = covariate-adjusted programs + calibrated sharing test + stability/extract ported with tests. (Core of that is done; remaining work is evaluation and the items below.)

## Done (committed on the branch)
- `2e88deb` `covariates=` (hard adjustment): per-cell-type `lm.fit` residuals over observed donors before centering; features chosen before adjustment.
- `2bfa994` covariate terms constant within a cell type are dropped for that cell type (one-level factors made `model.matrix` error; hit on COMBAT DC, which has no St_Georges donors).
- `913fd1c` `covariate_mode = "fixed"` (soft, SOFA-style): covariate columns as fixed flashier sample-side factors (`ebnm_normal` scale 10, OLS init), free programs fit jointly; `Y_used` = Y minus fitted covariate part.
- `d811eae` fallback ladder in `fit_celltype_programs`: point_normal, then no nullcheck; recorded in `fits[[ct]]$fallback` (flashier nullcheck hit a NaN ELBO on COMBAT DC, K=10 near residual df).
- `126753d`, `eb769e2` ported `principal_angles`, `sharing_spectrum`, `align_programs`, `match_programs`, `program_*`, `top_genes`, `assess_program_stability` (donor subsampling, not bootstrap). Program ids stay `<ct>_<k>`; no `fit$programs`/`fit$input`.
- `1194a12` calibrated pair test: statistic sum(cos^2), donor-permutation null, BH over callable pairs, `underpowered` flag (chance cosine `(sqrt(ka)+sqrt(kb))/sqrt(n) >= 1`; excluded from BH, never called). Null type-I 0.034 @0.05 (500 reps, n=40, k=3).
- `44f20a6` effect size `z_pair`, `excess_frac` (0 chance, 1 smaller subspace nested); needed because p saturates at 1/(n_perm+1).
- Tests: full suite passes (covariates 14, data-model 19, ebmf 16, extract-stability 7, mofacellular 13, principal-angles 4, sharing-calibration 13). Run without devtools: `R_LIBS_USER=/nonexistent /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin/Rscript -e 'library(stats);library(utils);for(f in list.files("R",full.names=TRUE)) source(f); library(testthat); source("tests/testthat/helper-toy.R"); test_dir("tests/testthat", env=globalenv(), reporter="summary")'`

## Results so far
- **COMBAT sharing** (K=10, scores space, 1999 perms; `benchmarks/biomni_replication/results/sharing/combat_repo_{none,Institute,Institute_Outcome}_{fit.rds,pairs.csv,pairs_v2.csv}`, scripts `07_adjusted_sharing.R/.sbatch`, `09_effect_size.R`):
  | adjustment | shared (q<.05) | underpowered | median excess_frac | median z | min z |
  |---|---|---|---|---|---|
  | none | 41 | 4 | 0.50 | 23.9 | 8.9 |
  | Institute | 40 | 5 | 0.49 | 23.6 | 6.6 |
  | Institute + Outcome | 39 | 6 | 0.36 | 17.6 | 5.3 |
  - All callable pairs still shared; every non-call is an underpowered pair, not a negative.
  - Institute changes nothing (removes program-batch correlation 0.2-0.5 -> ~0). Outcome removes ~a quarter of the effect.
  - Strongest, barely changed: cMono-ncMono, CD4-CD8, CD4-NK, CD8-NK, cMono-NK (within-lineage). Most reduced: pairs with DP and PB (B-DP 0.46 -> 0.14); may be adjustment cost (6 df on ~50 donors) not real loss.
  - Underpowered count rises because K is 10 on every cell type, not because sharing drops. DC used the no-nullcheck fallback (programs unchecked).
  - `Outcome` is coarse: blank for the 23 Sepsis donors (set to level `not_recorded` = effectively a Sepsis indicator in 07 script), mixes severity and disease type.
- **Simulation** (`08_sim_adjustment.R/.sbatch`, `summarize_sim_adj.R`, 20 reps, results in `results/sim_adj` (none/hard) and `results/sim_adj3` (none/hard/soft)): unadjusted test rejects 100% even with no sharing (batch alone); hard and soft both reject 0.15 / 0.00 with no sharing (rho 0 / 0.5) and 1.00 with sharing; recovery ~0.85 when the program is batch-confounded (rho 0.5). **Soft == hard** in every scenario. Decision rule was "soft must win on something hard cannot" -> it did not.

## Pending tasks (do in this order)
1. **Decide on the soft arm:** keep `covariate_mode = "fixed"` as an option, default stays `"residualize"`. Optionally run 100 reps for the 0.15 type-I (3/20, wide CI). Optionally compare against SOFA itself (`pip install biosofa` + muon in `cellprograms-py`, two cell types as views) as the item-6 comparator.
2. **More adjustment runs on COMBAT if wanted:** age bracket, sex, time since onset; user chose pure technical, then Institute + Outcome. Not run: Source/disease.
3. **Fix the K-vs-power problem:** K=10 everywhere makes small cell types uncallable; K selection or per-cell-type K (stability pruning) would help. DC fit is unchecked (fallback).
4. **Available-case representation** for concatenated programs (zero-fill leaks Institute; per-CT fits already observed-donor only).
5. **Default metrics:** CV ridge-logistic / per-program association tests with FDR plus dimension-matched nulls (`04_corrected_eval.R` has the floor and available-case kNN).
6. **Run `assess_program_stability` on COMBAT** (tool exists, never run on real data).
7. **Extensions:** reference mapping, case-control mixed models (dreamlet-style), cross-cohort matching with a null for |cosine|.
8. **Comparators not yet run:** MOFAcell, scITD, DIALOGUE, mc-ASTRA/MOFA-FLEX, SOFA.
9. **Not done from the replication:** `fit_joint_block` at K=10, GloScope/PILOT, cross-cohort matching; Stephenson/HLCA/OneK1K cp numbers vs Biomni (Biomni used K=25, we ran K=30).
10. **Commit** `benchmarks/`, `biomni/`, `HANDOFF.md`, `.gitignore` change (user's call; `.gitignore` already excludes `data/`, `results/fits/`, 26 GB of h5ad); push/PR the branch (user's call).

## Next action
Pick one of: (a) commit `benchmarks/` + `biomni/` on the branch (item 10), (b) `assess_program_stability` on `results/sharing/combat_repo_Institute_Outcome_fit.rds` (item 6, small K prune), (c) run SOFA as external comparator on the sim harness (item 1/8).

## State
- Branch `feat/covariates-and-sharing-tools`; uncommitted: `M .gitignore`; untracked `HANDOFF.md`, `benchmarks/`, `biomni/`
- Running: nothing of mine. Unrelated jobs of other users/projects in the queue (ejm-proj-seq, array, xbench-v2-final) are not mine to touch.
- SLURM: use `--partition=all` (QOS on `medium` caps 4 cpu/8G)

## Locked decisions
- Use `all` partition, not `medium` (observed, 2026-10-07)
- Repo defaults stay K=10, `point_laplace`, `var_type=1`; Biomni's K=30 understates cp (user approved plan, 2026-10-07)
- Do not port L1 shared priors, L2 init borrowing, EV-BIDIFAC, ICA/PVA (user, evidence in `biomni/REVIEW.md`)
- PMD not ported until residualization + calibrated null exist (`biomni/REVIEW_round2.md` §4); both exist now, so revisit
- Do NOT replace the repo's `fit_celltype_programs`; fork ported by adapter only
- First COMBAT adjustment: pure technical (Institute) (user, 2026-10-08); second: add `Outcome` (user, 2026-10-08)
- Hard adjustment (`"residualize"`) is the default; soft only if it wins on something hard cannot (user agreed framing, 2026-10-08)

## Dead ends — do not redo
- Pseudobulk from `raw.X` raw counts: does not reproduce Biomni. Correct input = sums of log-normalized `X`, floored to integers.
- Z-scoring the baseline representations: Biomni did not.
- Biomni's `max|Spearman|` with a 10-dim random floor: biased for ~266-dim reps; use dimension-matched floors.
- `cv_ridge` closed-form LOO without the intercept leverage term (fixed at `04_corrected_eval.R:57`).
- PA loadings-space for sharing among cell-type-specific genes: AUROC ~0.5-0.68.
- Bootstrap with replacement for stability (duplicated rows give junk low-rank factors); subsampling used.
- `ebnm_flat` as the prior for fixed covariate factors: no likelihood -> NA ELBO in flashier (use fixed-scale `ebnm_normal`).
- Using `p_pair` alone to compare COMBAT runs: saturates at 0.0005 for 41-44 pairs; use `excess_frac`/`z_pair`.
- Seeding inside `principal_angles` resets the RNG: do not pass a fixed `seed` in a loop that also generates data.

## Read first
- `biomni/REVIEW_round2.md` — all findings, replication table, PA vs PMD, port decisions
- `biomni/REVIEW.md` — round 1 (older; some claims revised)
- `R/principal-angles.R` (`.pa_engine`, `principal_angles`, `sharing_spectrum`), `R/ebmf.R` (`.cov_design`, `.residualize`, `.flash_fixed`, fallback ladder)
- `benchmarks/biomni_replication/README.md` — pseudobulk recipe, inferred choices
- `benchmarks/biomni_replication/results/sharing/*_pairs_v2.csv`, `results/sim_adj3/`

## Open items
- Saved COMBAT "repo_K10" fits in `results/fits/` were made with the FORK's `fit_celltype_programs` (ids `ct__F001`, have `programs`/`input`); the 07 fits use the repo function (ids `ct_k`). Compare 07 arms to each other, not to those.
- Replication residuals: OneK1K donor cap, clr pseudocount, `random` seed — only if exact Biomni parity matters
- A task "Review EGSV2 P vs NP findings F-003" appeared in the CellPrograms Biomni project (not from this session) — ignore unless the user asks
