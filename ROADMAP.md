# cellprograms — inventory and roadmap

**As of 2026-10-08.** Branch `feat/model-benchmark` @ `e524b21` (tag `prereg-v1`), PR #2 (`fix/review-defects`) open, benchmark pilot running.
**Rule for fit changes:** a modification enters `R/` only if it clears its pre-registered rule in `benchmarks/model_selection/PREREG.md`. Everything else keeps the current default.

## 1. Inventory

Status: **done** (merged or committed), **open PR**, **gated** (decided by the benchmark), **todo**, **decide** (needs your call), **rejected**.

### A. Package: already in place
| Feature | Where | Status |
|---|---|---|
| Per-cell-type EBMF fit, hard (residualize) and soft (fixed) covariates | `R/ebmf.R` | done (PR #1) |
| Fallback ladder, NaN-ELBO failure, per-attempt seed | `R/ebmf.R` | open PR #2 |
| `droplevels`, 10-df design guard, leverage warning | `R/ebmf.R` | open PR #2 |
| Scores-space sharing test, `z_pair`, `excess_frac`, underpowered flag, BH over callable pairs | `R/principal-angles.R` | done |
| Per-pair seeds, `n_perm` default scaled to number of pairs | `R/principal-angles.R` | open PR #2 |
| `strata=` blocked permutation null (D1) | `R/principal-angles.R` | committed on `feat/model-benchmark`, **no PR yet** |
| Stability refits keep `covariate_mode` / `flash_control` | `R/stability.R` | committed on `feat/model-benchmark`, **no PR yet** |
| Subsampling stability, `match_programs`, `align_programs`, `canonicalize_programs`, extract helpers, MOFAcellulaR export, simulator | `R/` | done |

### B. Fit-changing candidates (all gated on the benchmark)
| ID | Candidate | Arm | Rule it must clear |
|---|---|---|---|
| B1 | K by n: min(10, floor(n/10)) | `k_n` | dP on S3 ≥ +0.10 |
| B2 | K_max = 20 | `k_20` | dP on S3 ≥ +0.10 |
| B3 | K by stability (≥ 0.7 over 10 subsamples) | `k_stab` | dP on S3 ≥ +0.10 |
| B4 | Select genes after adjustment (D6) | `genes_after` | dR on S1+S2 ≥ +0.02 |
| B5 | `var_type = c(1,2)` | `vt_12` | dR on S1+S2 ≥ +0.02 |
| B6 | Cell-count weights via flashier `S` (D5) | `w_S` | dR on S1+S2 ≥ +0.02 |
| B7 | `strata=Institute` as the default null | `strata` | lowers F on S0 where base F > 0.05, P cost ≤ 0.05 |

All also need: CI excludes 0 on `semi`, same sign on `param`, F ≤ base + 0.03, dR ≥ −0.02, time ≤ 3× base.

### C. Package features not tied to the fit benchmark
| ID | Feature | Status | Note |
|---|---|---|---|
| C1 | **Coordinated-program detection across cell types**: program-level calls (which program in cell type a is coordinated with which in b, and across how many cell types) beyond measured donor covariates, returned as a program object / graph (D12) | **MAIN SCOPE** (you, 2026-10-08) | The pair-level test does not say which programs are coordinated. Needs its own pre-registered benchmark: program-level precision/recall on planted coordinated vs nuisance-driven programs, plus null calibration under a donor-level nuisance |
| C2 | Optimal (not greedy) program matching (D13a) | **decide** | Only worth it if it changes matches in practice; needs a small pre-registered check on stability refits (below) |
| C3 | Available-case downstream (D4) | todo, benchmarks only | Package already uses observed donors; zero-fill lives in `benchmarks/biomni_replication/02_metrics.R` |
| C4 | K ≤ n − q − 1 validation guard (D7c) | **decide** | Never binds at n ≥ 50, K ≤ 10, so it is a safety guard, not a performance change |
| C5 | Common gene panel across cell types (D6b) | deferred | Affects loadings space only; primary analysis is scores space |
| C6 | Docs for any new/changed argument (roxygen), tests per ported feature | todo, per port | |

### D. Analysis and benchmarks
| ID | Item | Status |
|---|---|---|
| D1 | Model-selection benchmark (harness, pre-registration, pilot, main run, decisions) | in progress |
| D2 | COMBAT rerun with final settings: per-pair seeds, per-cell-type K, `strata=Institute`, DC treated as sensitivity | todo, after M4 |
| D3 | `feat/sofa-comparator` (5 commits, tip 7efc410): SOFA harness, before/after recovery, ordinal severity | not pushed, no PR |
| D4 | Untracked results: `sim_adj4/`, SLURM logs under `benchmarks/biomni_replication/results/logs/` | todo (commit data worth keeping, ignore logs) |
| D5 | PMD comparison | not ported (principal angles ranked as well or better; PMD null rejects everything) |

### E. Rejected (do not revisit without new evidence)
L1 shared priors, L2 init borrowing, EV-BIDIFAC, ICA/PVA, replacing the repo's `fit_celltype_programs` (user decision, `biomni/REVIEW.md`); 7-level `Outcome` as primary adjustment; AUROC of `z_pair`; `p_pair` alone for comparing runs.

## 2. Roadmap

Dependencies: M0 → M1 → M2 → (M3) → M4 → M5. M6 and M7 are independent of M1–M5.

| # | Milestone | Deliverables | Exit criterion | Gate / depends on |
|---|---|---|---|---|
| **M0** | Land the fixes | Merge PR #2. Open PR for `feat/model-benchmark` (strata, stability fix, harness) | Both on `main`; full suite passes | You review the PRs |
| **M1** | Benchmark ready | Pilot done; `amplitudes.csv` and `DEVIATIONS.md` filled; per-dataset runtime measured | One amplitude per base × scenario chosen by the pre-registered rule; main-run cost confirmed (plan: about 45 min/dataset, 600 CPU-h, 6–10 h wall) | Pilot jobs finishing |
| **M2** | Main run and decisions | 900 datasets (450 per base) complete; `summarize.R` output `decisions.csv` and `decisions_base_calibration.csv` | Every arm has adopt / reject; `base` calibration and the confounder check reported; failed datasets rerun | M1 |
| **M3** | Confirmation (only if ≥ 2 arms adopted) | Combined winners vs `base` on seeds 10001+ | Same rules pass for the combination; otherwise only the best single arm is ported | M2 |
| **M4** | Port winners | New arguments or default flips in `R/`; one test each; docs; `HANDOFF.md` updated with locked decisions | Suite passes; benchmark arm reproduces from the new defaults; nothing ported that failed a rule | M2 (M3) |
| **M5** | COMBAT final analysis | Rerun with ported settings; per-pair seeds; DC sensitivity; compare to the old 0.50 / 0.49 / 0.36 / 0.43 `excess_frac` | Result table with CIs; every callable pair's call and every non-call stated with its reason | M4 |
| **M6** | Optional package features | C1, C2, C4 as you decide. C2 gets a pre-registered check: on stability refits, optimal vs greedy changes the matched set in ≥ 5% of programs, else drop it | Each accepted item has a spec, a test and a docs entry; rejected items recorded in `HANDOFF.md` | Your scope decision |
| **M7** | Housekeeping | PR or archive for `feat/sofa-comparator`; commit or ignore untracked results; C3 fix in benchmark metrics | No untracked benchmark outputs; SOFA branch either merged or documented as archived | Independent |

### What each milestone costs (rough)
- M0: minutes of your review time.
- M1: pilot is running now; picking amplitudes is a short script once it finishes.
- M2: 6–10 h wall time on partition `all`, depends on queue.
- M3: about 1 day wall including a new combined arm in `arms.R`.
- M4: about half a day per ported feature (code, test, docs).
- M5: about 10–30 min per COMBAT fit set on `all`, plus write-up.
- M6, M7: not estimated until scoped.

## 3. Open questions for you
1. ~~Is C1 in scope?~~ Yes, it is the main scope. Open: how are multi-cell-type programs defined (see the chat summary of PMD vs principal vectors)?
2. Should C4 (K ≤ n − q − 1) be added anyway as a guard, since it never changes results?
3. SOFA branch: open a PR, or archive it?
