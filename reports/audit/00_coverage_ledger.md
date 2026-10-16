# 00 — Coverage ledger (Phase 0)

Audit target: `mdmanurung/cellprograms` @ `413e290` (main), package 0.2.0. Baseline recorded 2026-10-17; working tree clean at audit start; no source modifications before this ledger. Environment: R 4.4.2 (conda), library `/workspace/.Rlib` — cellprograms 0.2.0, flashier 1.0.60 (willwerscheid@fd44811), flashier.utils 1.1.0.9000, ebnm 1.1.42, mashr 0.2.79, igraph 2.1.4. CI at `413e290`: R-CMD-check success, pages deploy success.

## Status vocabulary

- `source-reviewed`: code read in full this audit; runtime verification may still be pending.
- `runtime-launched`: dispatched to an independent verification workstream (see TASKS.md).
- `pending`: review not yet performed; covering phase listed.
- Risk: candidate severity per the audit prioritization (P0-P3) or `low`.

## 1. Core R implementation (17 files, all read)

| Component | Source | Scientific purpose | Tests | Reference | Status | Risk |
|---|---|---|---|---|---|---|
| Data contract | `R/data-model.R`: `as_cell_program_data`, `validate_cell_program_data` | Pseudobulk input model (per-cell-type donor x gene matrices + donor metadata) enforced before fitting | test-data-model.R (4) | Internal specification | source-reviewed; Inf/NaN contract gap suspected (CAND-03) | P2 |
| EBMF wrapper | `R/ebmf.R`: `fit_celltype_programs`, `.flash_fixed`, residualization, fallback ladder, union-panel fixed path | Per-cell-type EBMF via flashier; covariate adjustment (hard residualize or fixed mode); feature selection; seed control; fallback semantics (point_laplace -> point_normal -> no_nullcheck; non-finite ELBO = failure; all attempts re-seeded) | test-ebmf.R (6), test-covariates.R (4), test-sync-fitting.R (5), test-union-panel.R (4) | flashier (Willwerscheid et al.); ebnm (JSS 114(3)); Kotliar 2019 | source-reviewed; fixed-mode rank-deficiency crash suspected (CAND-02); parity vs direct flashier pending (child A) | P1 candidate |
| Preprocessing freeze | `R/preprocess-freeze.R` | Content-addressed (sha256) frozen preprocessing spec for exact replay | test-mr0-freeze.R (6) | Internal | source-reviewed; digest computed before self-reference; clean | low |
| Projection | `R/cell-projection.R`: `project_cells` | Fail-closed projection of frozen factorization onto new donors/cells; refuses legacy/fixed fits, undeclared feature_scale, design mismatch; sparse unsupported v1 (documented) | test-mr1-projection.R (7) | Internal + flashier-utils | source-reviewed; biological-validity checks pending (Phase 10) | P2 |
| Stability | `R/stability.R`: `assess_stability` | Donor-subsample refits; recovery frequency, loading/score similarity, sign consistency, top-gene overlap | test-extract-stability.R (3) | Internal; cf. cNMF stability practice | source-reviewed; failed-refit denominator issue (CAND-01); no `preprocessing_spec` replay (CAND-12) | P2 |
| Principal angles / matching | `R/principal-angles.R`: `principal_angles`, `sharing_spectrum`, `align_programs`, `match_programs`, `.pa_engine`, `.perm_within`, `.greedy_match`, `.pair_seed` | Subspace sharing tests between cell types (permutation null + random-subspace asymptotic reference + underpowered guard); factor matching across fits; deterministic per-pair seeds; n_perm floor warning vs BH level | test-principal-angles.R (4) | Internal geometry; random-subspace theory | source-reviewed; greedy-vs-"Hungarian-style" doc mismatch + NA guard (CAND-09); calibration pending (child D) | P2 |
| Program sharing | `R/program-sharing.R`: `assess_program_sharing`, `.sharing_matrices`, `.estimate_V`, `.match_programs_loadings`, `.degenerate_programs` | Shared/partial/private gene-loading classification via mashr EB covariance; CorShrink-style shrunk null-gene correlation; PSD scaling handled (`F_psd * scale_multiplier`); tiny-SE floor + broken-gene dropping documented | test-program-sharing.R (9), test-sharing-calibration.R (3) | mashr / Urbut 2019 | source-reviewed; gene-label permutation null under gene-gene covariance suspected anticonservative (CAND-07) | P1 candidate |
| CMP detection | `R/coordination.R`: `find_coordinated_programs`, `cmp_scores`, ambient scoring, permutation machinery, greedy graph | Cross-cell-type donor-level coordination: rank-normalized scores, strata residualization, Fisher-z T statistics, node aggregation A_v = sum_c max(0, max|T|-tau)^2, max-T permutation calibration (within-strata, per-cell-type joint row permutation), greedy graph, consensus scores | test-coordination.R (10) | Internal mathematical specification | source-reviewed; permutation design assessed sound for the complete null; link_threshold uncalibrated (CAND-04), sign-NA chains (CAND-05), singleton strata (CAND-06); calibration pending (child B) | P1 candidate |
| CMP interpretation | `R/cmp-interpret.R`: `prepare_collectri`, `annotate_cmp_tf`, `cmp_triage`, `cmp_nuisance_genes`, `cmp_context_association` | TF footprints (decoupler-style ULM + self-contained sign-flip gate), one-row triage table with nuisance flags, CollecTRI parsing (NFKB/AP1 merging), donor-variable association | test-cmp-interpret.R (4) | CollecTRI/OmniPath; decoupler ULM | source-reviewed; anticonservative competitive null documented in-code; VDJ regex misses joining genes (CAND-08); sign-NA default +1 propagates here | P3 |
| CMP annotation | `R/cmp-annotate.R`: `annotate_cmp_lr`, `annotate_cmp_themes`, `cmp_gene_overlap` | Ligand-receptor compatibility (degree-preserving resource permutation, FWER or pooled-FDR control, MOFAtalk-style donor check), theme ULM on |loading|, hypergeometric top-gene overlap classes | partial (test-cmp-interpret.R) | MOFAtalk-style check; LIANA-style resources | source-reviewed; dedicated tests sparse; causality language appropriately hedged in docs | P3 |
| Artifact genes | `R/artifact-genes.R` | Screening regexes: vdj, genotype (GST/HLA-DRB345/DQA/DQB, ERAP2, CCL3L1...), sex, IEG, platelet/RBC | test-artifact-genes.R (6) | Heuristic sets, documented as such | source-reviewed; CAND-08 | P3 |
| Extraction / canonicalization | `R/extract.R`: `canonicalize_programs`, top-gene helpers | Sign + RMS-scale canonicalization preserving reconstruction exactly (math verified by hand); PSD multipliers stored and applied downstream | test-ebmf.R, test-extract-stability.R | Internal | source-reviewed; runtime invariance + upstream flashier.utils parity pending (child A, CAND-11) | P2 |
| MOFAcellulaR export | `R/mofacellular.R` (88 ln) | Export pseudobulk as counts/coldata for MOFAcellulaR; requires complete cell_counts; TMM expects raw counts (documented) | test-mofacellular.R (2) | Ramirez Flores et al. 2023 (eLife 93161) | source-reviewed | low |
| Simulation | `R/simulation.R`: `simulate_pseudobulk` + scenarios (private/mixed/same_genes_independent/global) | Linear-Gaussian pseudobulk simulator with planted shared-activity vs shared-genes controls | test-simulation.R (1) | Internal | source-reviewed; LIMITATION: no count/NB/dropout/ambient process (CAND-10) | P2 (scope) |
| Plots | `R/plots.R` (8+ functions) | ggplot suite per analysis stage; graceful empty-plot degradation; membership alluvial with cross-run CMP label canonicalization | test-plots.R (4) | Internal | source-reviewed (remainder of file in saved tool output) | low |
| Utils / package | `R/utils.R`, `R/cellprograms-package.R` | Internal helpers (`.stopf`, `.warnf`, `%||%`), NSE globalVariables declaration | indirect | - | source-reviewed | low |

## 2. Benchmarks (`benchmarks/`)

| Benchmark | Size (tracked) | Purpose | Status | Risk |
|---|---|---|---|---|
| batch_adjustment | 24K | ComBat / ComBat-seq / residualization comparison | pending (Phases 4, 12): comparator fairness when simulation matches one comparator's assumptions | P2 |
| biomni_replication | 8.8M | CELLxGENE Census pseudobulk replication of external pipeline (raw-count vs transformed representation question) | pending (Phases 2, 13) | P2 |
| c2_matching | 9.5K | Cross-cohort matching artifacts | pending (Phase 12) | P3 |
| cmp_bench | 259K | CMP detection benchmark development artifacts | pending (Phase 12) | P2 |
| cmp_detection | 535K | cmp-bench-v1.2 preregistered held-out confirmation (pooled membership F1 0.571 -> 0.705 at cap 30; paired diff 0.127, CI 0.113-0.142) | pending (Phase 12): effect decomposition + held-out reuse question | P1 |
| coordination | 58K | Coordination benchmark artifacts | pending (Phase 12) | P2 |
| divas_comparison | 70K | DIVAS reanalysis comparison | pending (Phase 13) | P2 |
| model_selection | 72K | Rank-selection evidence (K=10 vs K>=30 operating points) | pending (Phase 3D) | P2 |
| mofacellular_comparison | 58K | MOFAcellulaR comparison | pending (Phases 11, 13) | P2 |

## 3. Tests (`tests/testthat/`)

20 test files, 100 `test_that` blocks (plus 2 helpers). Classification into API-contract / mathematical-invariant / numerical-reference / statistical-calibration / biological-negative-control / integration categories: pending (Phase 16). Known coverage from source reading: end-to-end fit + canonicalization reconstruction invariant (test-ebmf.R), degenerate K=0 sentinel, freeze/replay, projection fail-closed, sync parity fixtures. Local suite blocker: testthat not installed in sandbox R libraries (CI-only); installation dispatched.

## 4. Vignettes, docs, python, analyses

- `vignettes/`: quickstart + 4 articles (real_data_case_study, immune_aging_cross_cohort, divas_case_study, mofacellular_comparison); reconciled with `docs/` in the prior session (commit `413e290`); content-level scientific audit pending (Phases 13, 17).
- `docs/`: committed pkgdown 2.2.1 site (GitHub Pages, legacy build); URL-checked (6/6 correct).
- `python/`: 25 scripts (COMBAT baselines/metrics/panels/regimes, MOFAflex runs, patpy comparison, Stephenson runs, tuning search, vartype2, winner characterization) — R/Python boundary audit pending (Phase 15).
- `analyses/immune_aging/`: cross-cohort aging analysis (OneK1K 981, AIDA 625, CIMA 421, ABF300 166, KCL 188, UCSF 99, SoundLife 96, SPAC 33 = 2,609 donors) — pending (Phase 13).
- Reproducibility bundle: `mdmanurung/cellprograms-reproducibility` @ `ec9bf2d` (505 files), pins `cellprograms@413e290`.

## 5. Reports and release evidence

| Report | Lines | Content | Reconciliation status |
|---|---|---|---|
| `reports/ROADMAP.md` | 147 | Inventory + roadmap; Stage A closed 2026-10-16 | pending (Phase 11) |
| `reports/STAGE_B_ROADMAP.md` | 95 | MR-0 freeze + MR-1 projection IMPLEMENTED; MR-2 cNMF, MR-3 metacells, MR-4 mcRigor/anchored STAGED | pending (Phase 11); staged work must not be presented as complete |
| `reports/SYNC_ROADMAP.md` | 112 | Stage A sync receipt: dependency direction, supported version pair, parity evidence, known deviations, rollback | pending (Phase 14) |
| `reports/plan_review.md` | 87 | Plan review (statistical soundness, benchmark gaps, feasibility) | pending (Phase 12) |
| `reports/final_report.md` | 125 | Attempt history, diagnostics, ablations, leakage audit | pending (Phases 12, 21) |
| `reports/CRAN_PREP_CHECKLIST.md` | 46 | Submission checklist | pending (Phase 18) |
| `benchmarks/cmp_detection/REPORT_v2.md` | 180 | v2 report incl. held-out results, deviations, limitations | pending (Phase 12) |
| `benchmarks/cmp_detection/PREREG_cmp-bench-v1.2.md` | 35 | Pre-registration: frozen candidate, comparators, criteria, prohibitions | pending (Phase 12): chronology audit |

## 6. Candidate findings register

See `TASKS.md` (CAND-01 … CAND-12). None reclassified yet; every candidate requires a reproducer plus independently established expected behavior before CONFIRMED status. No P0 candidates identified at source level so far; the two P1 candidates (CAND-04, CAND-07) concern statistical calibration of downstream inference, not core factorization mathematics.
