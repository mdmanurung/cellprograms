# Handoff — cellprograms: cross-cell-type sharing, model-selection benchmark, coordination methods

**Updated:** 2026-10-10 · **Checkout:** `feat/coordination` = local `main` = `origin/main` · **Status:** M0+M1 and all five benchmark fixes complete; **M2 main run finished and summarized (2026-10-10): no arm adopted under the registered rules** (see `benchmarks/model_selection/RESULTS.md`); the `strata` verdict waits on a user decision. `feat/sofa-comparator` and `feat/model-benchmark` remain unpushed.

> Note on this file: an earlier version's "Next action" (push `fix/review-defects`) and the `strata=` blocked-null work are now **complete and on `main`** — see "Recently completed" below. This file was refreshed rather than creating a near-twin `handoff.md`, which would collide on case-insensitive filesystems.

## Revised gene-sharing plan — saved 2026-10-09

The user requested a local revised implementation plan after a source and statistical review. The reviewed progress tracker is [docs/plans/2026-10-09-selective-coupled-ebmf-gene-sharing.md](docs/plans/2026-10-09-selective-coupled-ebmf-gene-sharing.md). It defines a **separate gene-loading recurrence workstream**, preserving independent `fit_celltype_programs()`, the coordination estimand, and the existing `prereg-v1` benchmark plus recorded deviations.

- **Next step for this plan:** GS-00, a specification-only model and selection contract. Exact penalties/precision, sharing thresholds, eligibility, rank/ambiguity rules, recruitment, and reporting must be resolved before the optimizer. GS-05 and later are explicitly blocked on their prerequisites.
- **First comparison:** include signed EBMF consensus **plus fixed-template residual recruitment** (B1-R); test weak true members and genuine nonmembers. Joint template updating must demonstrate incremental benefit over that comparator before flexible or empirical Bayes extensions.
- **Compatibility work tracked:** MOFA gene-row alignment, stability success/failure denominators, optimal benchmark assignment, coherent gene-loading clusters, and dependency-enabled package checks. These repairs were reviewed, not implemented in this session.
- **Review evidence:** dependency-free checks reproduced the adapter/matching/clustering/stability issues; both coordination smoke scripts passed. The R4_51 package-test attempt had 15 skips and one dependency error because `flashier`/`ebnm` were unavailable. No complete passing package suite, hierarchical implementation, calibration, or HPC benchmark result is claimed.
- **Authorization/state:** this request authorized documentation only. The plan and this note were written; package/benchmark source and the recorded active main run were not changed or resubmitted. Queue counts and runtime state below are historical observations, not refreshed by this documentation update. Keep runtime-source files immutable while existing jobs remain active; future implementation needs its own action request and an isolated checkout.

The reviewed hash `584c50236c78730f4b33fe2d82bfd211b71a5629` is the commit ID; its Git tree is `f975b9ec1e4b074db762359f261db9628e6397bd`. The saved plan is **not yet preregistered**, and no implementation step is complete. Existing operational handoff content below is retained.

## Goal
Decide which fit modifications earn their complexity (model-selection benchmark), and land a defensible cross-cell-type sharing/coordinated-program capability (coordination benchmark, roadmap item C1). Both are pre-registered; nothing enters `R/` unless it clears its registered rule.

## Recently completed (all on `main`)
- PR #1 (covariates, calibrated sharing test, extract/stability) and PR #2 (review defects) merged.
- `strata=` blocked permutation null in `principal_angles`/`sharing_spectrum` (`.pa_engine`, `R/principal-angles.R`).
- Model-selection harness frozen at tag `prereg-v1` (`benchmarks/model_selection/`): arms `k_n`, `k_20`, `k_stab`, `genes_after`, `vt_12`, `w_S`, `strata`, diagnostics `none`/`none_strata`.
- Amplitudes chosen by pilot (`amplitudes.csv`, `DEVIATIONS.md` D1-D2). Pilot runs only `base`; the runtime pilot (all arms) measured per-arm fit time.
- Coordination simulator `sim_coord.R` + self-check (all pass). Candidate methods `methods.R` committed but `node_cor` broken (defect D below).

## Next action — decide the `strata` verdict, then pick the metaprogram direction
M2 is complete: 900 datasets, 8 arms, no failed tasks, `summarize.R` run, `decisions*.csv` and `RESULTS.md` committed. **No arm is adopted** as registered; M3 is not triggered. `strata` clears every rule except the D4 `frac_untested` guard, which fails every `param` arm including `base` (12-27% untested). **Open decision (user):** keep `strata` rejected as registered, or report it as adopted under a post-hoc `base`-relative guard, logged as a deviation made after seeing results. The M5 COMBAT rerun uses `strata=Institute` only if `strata` is adopted. The confounder check passes on `param` (`none` F = 1.0) and fails on `semi` (`none` F = 0.09).

Metaprogram direction: nothing is chosen or implemented in `R/`. Options and the staged plan are in `docs/metaprogram_methods_review.md`; the coordination prereg still needs its four blockers fixed (method list, replicates and paired inference, node-to-truth matching rule, pinned EBMF config). Another session saved a gene-sharing plan above and has a separate clone `.c1-work/` (hidden locally via `.git/info/exclude`) with its own unpushed C1 commits; do not delete it.

## Current state
- Tests: latest Claude run reported 155 expectations, 0 failures after the namespace-mocking fix (`09f1fcc`). The 2026-10-09 progress check validated outputs and logs; it did not rerun package tests.
- `.gitignore` (on main) now excludes regenerable benchmark bulk (per-dataset `rep*.csv`, simulated matrices, fits, SLURM logs) while keeping scripts, pre-regs, and summary tables; `summarize.R`'s `decisions*.csv` at the model_selection results root stays trackable.
- Running: nothing of mine (M2 finished 2026-10-10). Use SLURM `--partition=all`.
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
**Status 2026-10-09, pushed to `origin/main`:** all 5 commits done (`613bc7c`, `1272736`, `1ae3908`, `73bbdd4`, plus `49dd869`, the `node_cor` fix; DEVIATIONS D3-D6). The `node_cor` fix: `.perm_rows` now keeps the original row names, so name-indexed callers see permuted values; with it node_cor recovers R100/R50/R20 and all 4 methods return 0 clusters on the null. `test_methods.R` now asserts this and fails on the old code. Smoke re-run: `decisions.csv` written, absent arms degrade to NA / not adopted. Full suite: 155 expectations, 0 failures after `09f1fcc` fixed the stability test's namespace mock. `584c502` added the array-chunk offset (D7) and was pushed before launch.

| # | Commit | Files | Nature |
|---|---|---|---|
| 1 | `fix(benchmark): record untested sharing tests instead of scoring them negative` | `arms.R:88-100` | additive `tested` + reason columns |
| 2 | `fix(benchmark): gate adoption on sharing-test coverage` | `summarize.R`, `DEVIATIONS.md` | additive guardrail, logged |
| 3 | `fix(principal-angles): clear error when no cell type retains a program` | `R/principal-angles.R:278`, `test-principal-angles.R` | message-only |
| 4 | `fix(coordination): make the node_cor permutation donor-level (positional)` | `methods.R`, `test_methods.R` | draft/unregistered |
| 5 | `docs(benchmark): record the 8-arm main run` | `DEVIATIONS.md` | deviation log |

Validate after each: full `testthat` suite; `test_methods.R` end-to-end; re-run the M2 smoke to confirm `decisions.csv` still produced and absent arms degrade to `NA`/not-adopted (not an error). Commit and push to `main` only once green.

## Main run (complete)
- Submitted with `ARMS=base,strata,none_strata,none,k_n,k_20,vt_12,genes_after` (drops `k_stab`, `w_S`).
- 900 datasets (450 per base), submitted in 10 array chunks (D7), using `tasks_param.txt` / `tasks_semi.txt`. Outputs land in gitignored `benchmarks/model_selection/results/<base>/<scenario>/`.
- Done: results, `base` calibration with Clopper-Pearson CIs and the confounder check are in `benchmarks/model_selection/RESULTS.md`; per-arm verdicts in `results/decisions.csv`.

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
