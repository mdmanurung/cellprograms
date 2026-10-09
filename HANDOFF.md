# Handoff — cellprograms: cross-cell-type sharing, model-selection benchmark, coordination methods

**Date:** 2026-10-08 · **Branch:** `main` @ `bdb8c69` (all branches pushed except `feat/sofa-comparator`) · **Status:** M0+M1 done (PR #1,#2 merged; harness, amplitudes, DEVIATIONS on main); M2 main run not yet started; three benchmark defects found by an end-to-end smoke test and diagnosed below, fix plan ready but not yet implemented.

> Note on this file: an earlier version's "Next action" (push `fix/review-defects`) and the `strata=` blocked-null work are now **complete and on `main`** — see "Recently completed" below. This file was refreshed rather than creating a near-twin `handoff.md`, which would collide on case-insensitive filesystems.

## Goal
Decide which fit modifications earn their complexity (model-selection benchmark), and land a defensible cross-cell-type sharing/coordinated-program capability (coordination benchmark, roadmap item C1). Both are pre-registered; nothing enters `R/` unless it clears its registered rule.

## Recently completed (all on `main`)
- PR #1 (covariates, calibrated sharing test, extract/stability) and PR #2 (review defects) merged.
- `strata=` blocked permutation null in `principal_angles`/`sharing_spectrum` (`.pa_engine`, `R/principal-angles.R`).
- Model-selection harness frozen at tag `prereg-v1` (`benchmarks/model_selection/`): arms `k_n`, `k_20`, `k_stab`, `genes_after`, `vt_12`, `w_S`, `strata`, diagnostics `none`/`none_strata`.
- Amplitudes chosen by pilot (`amplitudes.csv`, `DEVIATIONS.md` D1-D2). Pilot runs only `base`; the runtime pilot (all arms) measured per-arm fit time.
- Coordination simulator `sim_coord.R` + self-check (all pass). Candidate methods `methods.R` committed but `node_cor` broken (defect D below).

## Next action — implement the fix plan (5 commits, then push `main`)
See "Fix plan" section. Execute in order; validate after each. Then run the M2 main run (see "Main run" section).

## Current state
- Tests: full suite green (covariates 14, review-fixes 13, sharing-calibration 13, ebmf 16, data-model 19, extract-stability 7, mofacellular 13, principal-angles 4, simulation 4). Run command in the commit history / this file's earlier version.
- `.gitignore` (on main) now excludes regenerable benchmark bulk (per-dataset `rep*.csv`, simulated matrices, fits, SLURM logs) while keeping scripts, pre-regs, and summary tables; `summarize.R`'s `decisions*.csv` at the model_selection results root stays trackable.
- Running: nothing of mine. Use SLURM `--partition=all`.
- `feat/sofa-comparator` (5 commits, tip `7efc410`): SOFA comparator, before/after recovery, ordinal severity — **unpushed, no PR** (roadmap D3 / M7 open question).

## M2 smoke-test findings (end-to-end path verified)
Ran 6 datasets (param/semi × S0/S1/S7) with `MS_SMOKE=1`. Confirmed `summarize.R` runs end-to-end and writes `decisions.csv` + `decisions_base_calibration.csv` (it had never been exercised; it simply needed S0/S7, which the pilot lacks). The multi-day run is safe to launch. This surfaced three defects (below) plus the `k_stab`/`w_S` time-guard fact.

### Time guard already disqualifies two arms (prereg: median time ≤ 3× `base`)
From the full-scale runtime pilot (all 10 arms, 1 dataset/base): `k_stab` **9.73×**, `w_S` **4.51×**; `k_20`/`vt_12` ~1.25×, `k_n`/`genes_after` <1×. `summarize.R` sets `guards_pass=FALSE` above 3×, so `k_stab` and `w_S` cannot be adopted regardless of effect. **Decision (user): drop `k_stab` and `w_S` from the main run** — done via the `ARMS` submit-time env override in `run.sbatch`/`run_rep.R` (no frozen-file change), logged in `DEVIATIONS.md`. 8 arms → ~324 CPU-h (~1 day at `%14`) vs ~1006 CPU-h (~3 days) for all 10.

### Defect A (critical) — untested sharing tests scored as negative
`arms.R:88-100`: `sharing_spectrum` is wrapped in `tryCatch`; on error `ss=NULL`, so `get("shared", FALSE)` returns `FALSE` for every pair. A degenerate arm scores `F=0` — indistinguishable from perfect calibration. `F` gates guardrail 1 for all arms and is the primary endpoint for `strata`. **Fix:** add `pair_rows$tested` (`FALSE` when `ss` is NULL or the pair is absent) + a failure-reason column; untested is never a tested negative.

### Defect B — no guardrail on test coverage
**Fix (user chose option (a)):** keep the preregistered `F` definition unchanged; add `frac_untested` per `(arm, base)` as a guardrail `≤ 0.05` + a column in `decisions.csv`. Log the additive guardrail in `DEVIATIONS.md`. A broken arm is rejected on coverage (and independently on effect/time), not rewarded with a fake `F=0`.

### Defect C — opaque `n < m` error
`R/principal-angles.R:278-279`: if no cell type retains a program, `combn(character(0), 2L)` throws `n < m`. **Fix:** explicit `length(cts) < 2` guard with a named message; message-only, no numeric effect; add a case to `tests/testthat/test-principal-angles.R`.

### Defect D — `node_cor` permutation is a no-op (coordination/C1)
`benchmarks/coordination/methods.R`: `stat_fun` re-indexes rows **by donor name** after permuting, restoring the original donor→value mapping (observed max-T == null max-T == 71.44685 in all 199 draws → every `fwer_p=1`, zero clusters; smoke test aborted). Scope confirmed confined to `node_cor`; `.pa_engine`/`sharing_spectrum`, `method_subspace_sum`, `method_gene_match` all align positionally and are correct. **Fix:** make `.maxT_edges`'s `perm_fun` pair-specific; pre-align each `node_cor` pair to positional matrices and permute positionally within strata; adapt `gene_match` to the new signature. `subspace_sum` unaffected.

### Defect E — `test_methods.R` asserted nothing (why D reached main)
The smoke test only prints; its sole failure was an incidental `aggregate` error. **Fix:** make `show()` robust to zero clusters and add a hard assertion that all 4 methods return **0 clusters on the null** dataset — a real regression guard against exactly this no-op-permutation bug.

### Not changing: `arm_k_stab`'s 0.7 stability threshold
At the preregistered amplitudes (`base` power ≈ 0.5, weak by design) `recovery_freq` is < 0.7 for every program, so `arm_k_stab` legitimately prunes to K=0. Adding a "keep ≥1 program" floor would rescue an arm the prereg defined at 0.7 and change B3's meaning — a deviation not recommended. Correct response is A+B: let it degenerate, report uncovered, reject with stated reason. With `k_stab` dropped from the main run this is moot for M2 but the coverage guardrail still protects other arms.

## Fix plan (5 commits, validated after each)
**Status 2026-10-09, on `feat/coordination`, not pushed:** commits 1, 2, 3, 5 done (`613bc7c`, `1272736`, `1ae3908`, `73bbdd4`; DEVIATIONS D3-D6). Commit 4 (`node_cor` positional permutation + real assertions) still open; amend per the Biomni review (one shared permutation draw per iteration, positional strata, pinned seed). Smoke re-run: `decisions.csv` written, absent arms degrade to NA / not adopted. Full suite: 1 failure, `test-extract-stability.R` "stability refits keep covariate_mode and flash_control", also fails on the untouched tree (pre-existing; the "suite green" line below is out of date).

| # | Commit | Files | Nature |
|---|---|---|---|
| 1 | `fix(benchmark): record untested sharing tests instead of scoring them negative` | `arms.R:88-100` | additive `tested` + reason columns |
| 2 | `fix(benchmark): gate adoption on sharing-test coverage` | `summarize.R`, `DEVIATIONS.md` | additive guardrail, logged |
| 3 | `fix(principal-angles): clear error when no cell type retains a program` | `R/principal-angles.R:278`, `test-principal-angles.R` | message-only |
| 4 | `fix(coordination): make the node_cor permutation donor-level (positional)` | `methods.R`, `test_methods.R` | draft/unregistered |
| 5 | `docs(benchmark): record the 8-arm main run` | `DEVIATIONS.md` | deviation log |

Validate after each: full `testthat` suite; `test_methods.R` end-to-end; re-run the M2 smoke to confirm `decisions.csv` still produced and absent arms degrade to `NA`/not-adopted (not an error). Commit and push to `main` only once green.

## Main run (after the fix plan)
- Submit with `ARMS=base,strata,none_strata,none,k_n,k_20,vt_12,genes_after` (drops `k_stab`, `w_S`).
- 900 datasets (450 per base). `sbatch --array` per `tasks_param.txt` / `tasks_semi.txt`. Outputs land in gitignored `results/<base>/<scenario>/`.
- After completion: `summarize.R <results_root>` → `decisions.csv`; report every arm adopt/reject, `base` calibration + Clopper-Pearson CI, and the confounder check (`none` must show inflated F on S0 or the `strata` result is uninformative).

## Locked decisions
- **Coverage guardrail over redefining F** (user, 2026-10-08): keep preregistered `F`; add `frac_untested ≤ 0.05`.
- **Main run drops `k_stab` and `w_S`** (user, 2026-10-08): both breach the 3× time guard; `k_stab` also degenerates to K=0.
- Hard adjustment (`covariate_mode="residualize"`) default; soft an option only.
- Do not port L1/L2 priors, EV-BIDIFAC, ICA/PVA; do NOT replace repo `fit_celltype_programs`.
- Repo defaults K=10, `point_laplace`, `var_type=1`. Partition `all`, not `medium`.

## Dead ends — do not redo
- `summarize.R` "fails" on pilot data with `no rows to aggregate` — it needs S0/S7, which the pilot lacks. Not a bug.
- `p_pair` alone to compare runs (saturates); use `excess_frac`/`z_pair`.
- `ebnm_flat` for fixed covariate factors (NA ELBO); use fixed-scale `ebnm_normal`.
- AUROC of `z_pair` as evidence (1.00 even unadjusted); report FPR/power.
- Bootstrap with replacement for stability; fixed seed inside a data-generating loop.
- Adding a floor to `arm_k_stab` to rescue it — rejected; let coverage guardrail reject it.

## Read first
- `benchmarks/model_selection/PREREG.md` (frozen at `prereg-v1`) and `DEVIATIONS.md`.
- `benchmarks/model_selection/arms.R:57` (`arm_k_stab`), `arms.R:88-100` (Defect A), `summarize.R:38-42,69` (Defect B).
- `R/principal-angles.R` (`.pa_engine`, `sharing_spectrum`), `benchmarks/coordination/methods.R` (Defect D).

## Open items
- M6/C1 (main scope): multi-cell-type program object / program graph (D12), after `node_cor` is fixed. Needs its own pre-registered benchmark (see `ROADMAP.md`).
- `feat/sofa-comparator`: open a PR or archive it (roadmap D3 / M7).
- Not fixed from the Biomni review: cell-count weighting via flashier `S` (D5); common gene panel / post-adjustment gene selection (D6); K cap ≤ n−q−1 (D7c); optimal (not greedy) matching (D13a); available-case downstream (D4).
