# Audit task ledger — adversarial scientific audit of cellprograms

Persistent ledger. If execution is interrupted, resume from this file. Updated at each milestone.

## Baseline (frozen)

- Repo: `mdmanurung/cellprograms` @ `413e290` (main), v0.2.0; working tree clean at audit start.
- `flashier-utils` @ `271f583` (unchanged). Reproducibility bundle: `mdmanurung/cellprograms-reproducibility` @ `ec9bf2d`.
- R 4.4.2 (conda) on worker-1; library `/workspace/.Rlib`: cellprograms 0.2.0, flashier 1.0.60 (willwerscheid@fd44811), flashier.utils 1.1.0.9000, ebnm 1.1.42, mashr 0.2.79, igraph 2.1.4.
- CI at `413e290`: R-CMD-check success, pages deploy success.
- Rules in force: no source modifications before Phases 0-1 sign-off; every confirmed defect follows demonstrator -> independent expected -> regression test -> minimal fix -> full suite -> document changed science; benchmark acceptance criteria never changed post hoc.

## Status vocabulary

pending | in-progress | done | blocked. Finding classes: CONFIRMED / HYPOTHESIS / LIMITATION / VERIFIED-CORRECT / FALSE-POSITIVE.

## Phase tracker

| Phase | Scope | Status | Evidence |
|---|---|---|---|
| 0 | Reconnaissance + coverage ledger | done | `reports/audit/00_coverage_ledger.md` |
| 1 | Specification-to-implementation matrix | in-progress | `reports/audit/01_spec_matrix.md` |
| 2+4 | Pseudobulk/data validity + confounding | in-progress (child C + inline) | `shared/audit/C_stab/` |
| 3+5+6 | EBMF numerics + canonicalization + stability | in-progress (children A, C) | `shared/audit/A_ebmf/`, `C_stab/` |
| 7+8 | CMP inference + program sharing | in-progress (children B, D) | `shared/audit/B_cmp/`, `D_sharing/` |
| 9+10 | Interpretation + projection | pending | inline |
| 11+12+14 | Benchmarks + real data + sync | pending | inline |
| 15-18 | Software/tests/UX/publication | pending | inline |
| 19-21 | Panel + remediation + deliverables | pending | inline |

## Verification workstreams (background)

| ID | Scope | Status | Evidence dir |
|---|---|---|---|
| A | EBMF numerical validation: reconstruction invariance, parity vs direct flashier, seed reproducibility, zero/near-zero factors, var_type behavior, fixed-mode rank deficiency (CAND-02), posterior-SD transforms, ELBO/fallback semantics | launched | `shared/audit/A_ebmf/` |
| B | CMP calibration: complete-null FWER, partial-alternative strong control + membership precision (CAND-04), sign-NA chains (CAND-05), singleton strata (CAND-06), anticorrelation | launched | `shared/audit/B_cmp/` |
| C | Stability + data contracts: failed-refit denominator (CAND-01), Inf contract (CAND-03), edge cases, preprocessing_spec replay gap (CAND-12), metric response to signal/noise | launched | `shared/audit/C_stab/` |
| D | Sharing matching: gene-block permutation null (CAND-07), planted-truth classification, greedy-vs-optimal assignment (CAND-09), sign handling, .estimate_V PSD | launched | `shared/audit/D_sharing/` |

## Findings registry

| ID | Component | Claim | Class | Priority | Status |
|---|---|---|---|---|---|
| CAND-01 | R/stability.R | `recovery_freq = hits/n_boot`; failed refits silently count as non-recoveries; no `n_failed` exposed | code-confirmed behavior; class pending runtime demo | P2 | child C |
| CAND-02 | R/ebmf.R `.flash_fixed` | `solve(crossprod(Xc))` crashes on rank-deficient covariate design (residualize path has a rank guard; fixed path does not) | HYPOTHESIS | P2 | child A |
| CAND-03 | R/data-model.R | `any(is.na(m))` misses `Inf`; non-finite values pass the contract | HYPOTHESIS | P2 | child C |
| CAND-04 | R/coordination.R | `link_threshold=3` edge admission uncalibrated; membership precision not FWER-controlled under partial alternatives | HYPOTHESIS | P1 | child B |
| CAND-05 | R/coordination.R | chain-joined members can have `sign=NA`; `cmp_scores` defaults NA->+1 (arbitrary consensus orientation); propagates to annotate/interpret layers | code-confirmed behavior; impact pending | P2 | child B |
| CAND-06 | R/coordination.R | singleton strata -> leverage-1 residuals ~0 in `.cmp_prep` without warning | HYPOTHESIS | P3 | child B |
| CAND-07 | R/program-sharing.R | gene-label permutation null ignores gene-gene covariance; potentially anticonservative under correlated gene blocks | HYPOTHESIS | P1 | child D |
| CAND-08 | R/artifact-genes.R, R/cmp-interpret.R | vdj regex misses IGHJ/IGKJ/IGLJ joining genes (and IGLL5/IGLL1, TR*J); sets documented as screening heuristics | code-confirmed; low severity (documented heuristic) | P3 | open |
| CAND-09 | R/principal-angles.R | `match_programs` uses greedy matching; docstring says "Hungarian-style" (misleading); `.greedy_match` errors on NA input | code-confirmed doc mismatch; NA crash pending | P3 | child D |
| CAND-10 | R/simulation.R | linear-Gaussian only; no count/NB/dropout process | LIMITATION | P2 (scope) | documented |
| CAND-11 | R/extract.R | upstream flashier.utils canonicalization parity claimed "verified to 0 difference" but independently unverified | HYPOTHESIS | P3 | child A |
| CAND-12 | R/stability.R | refit replays feature sets + covariates + seed but not `preprocessing_spec` (provenance gap) | code-confirmed | P3 | child C |

## Decision log

- 2026-10-17: Baseline frozen at `413e290`; audit evidence collected strictly read-only.
- 2026-10-17: Verification children A-D dispatched against the pinned environment; repo untouched during evidence collection.
- 2026-10-17: Comparator citations verified (see `01_spec_matrix.md`): GeneNMF (Andreatta & Carmona, CRAN v0.9.2), MOFAcellulaR / multicellular factor analysis (Ramirez Flores et al. 2023, eLife 93161), DIVAS (Sun et al. 2026, Bioinformatics btag615; Prothero et al. 2022, arXiv:2212.00703), mcRigor (Liu & Li 2025, Nat Commun s41467-025-63626-5), Spectra (Kunes et al. 2023, Nat Biotechnol s41587-023-01940-3).
- 2026-10-17: testthat absent from sandbox R libraries (CI-only dependency); installation dispatched so the full suite can run locally during remediation.
