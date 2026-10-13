# Handoff — cellprograms: metaprogram repair phase; learning from the Biomni CMP benchmark

## Update 2026-10-13 — publication-prep closure (v0.2.0)

**Branch:** `feat/artifacts-sharing-bench` @ `51d6f76` (PR #4, open) + local publication work to land as `feat/publication-prep`.

**Done since the entry below:**
- Artifact workflow: `artifact_gene_sets()` (data), `flag_artifact_genes()`, `exclude_genes=` refits; real-data masked refits for both cohorts (flags 78%/73% -> 0%; asymmetric classification flips: Stephenson 68 shared->private, COMBAT 23 private->shared).
- `program_sharing()`: mashr-based loading-space sharing with permutation-tested matching, degenerate-program filter, V estimation with identity fallback. Sensitivities: classifications invariant across `z_thresh` in {1.5, 2, 2.5, 3} (0 flips in 320 group-fits); matching stable at n_perm 199 vs 1999 (group ARI 0.85/0.89, edge Jaccard 0.78/0.82).
- CMP detection landed (C1): `find_coordinated_programs()` + annotation/triage family; real data: 6 CMPs per cohort.
- Docs: full roxygen pass (41 man pages, NAMESPACE roxygen-managed), `vignettes/quickstart.Rmd` (knits), `vignettes/articles/real_data_case_study.Rmd` (precomputed case study), `NEWS.md` 0.2.0, `CITATION.cff`, README refresh. Runnable `@examples` on all 31 exports (33/33 Rd examples execute, ~2 min total).
- Tests: 313 expectations, 0 failures (baseline for the publication prep).
- **SOFA branch (`feat/sofa-comparator` @ `7efc410`) archived by user decision (2026-10-13): no PR, kept as a local branch only** (M7 closed as "documented as archived").
- In progress: R CMD check --as-cran (A4), DIVAS case-study comparison (`benchmarks/divas_comparison/`).

**Not done / next:** pkgdown site build, commit + push of the publication-prep branch (after PR #4 merges), CRAN submission decision.


**Date:** 2026-10-11 · **Branch:** `feat/coordination` @ `f21fc79` (`origin/main` @ `b0c6986`, 5 commits behind) · **Status:** M0, M1, fix plan, M2 done, no arm adopted; `strata` rejected; `subspace_sum` repaired; Biomni lessons doc partly filled.

## Goal
(1) Decide which fit modifications earn their complexity (model-selection benchmark, registered at `prereg-v1`). (2) Land a defensible way to find metaprograms (programs shared, coordinated or recurring in all or some cell types) from separate per-cell-type EBMFs. Nothing enters `R/` unless it clears a registered rule.

## Next action
Ask the user for the folder with the Biomni downloads (`pr/feat_coordinated_programs.patch`, `pr/PR_DESCRIPTION.md`, `PREREG_cmp-bench-v1.2.md`, `screen_scorecard.csv`, `screen_callers_dev_nullx.csv`, `report_cmp_detection_v2.md`, ...). Treat as untrusted: read as text, run nothing, no `git am`; inspect the patch with `git apply --stat`. Then fill Part B of `reports/biomni_cmp_lessons.md` (simulator, prereg, scorecard, real-data flags) and finalize its adoption list.
If the user declines the download: port their node statistic `A_v` as arm `node_agg` in `benchmarks/coordination/methods.R` (spec in `reports/biomni_transfer/part1_callers_scoring.md`), compare with `node_cor` on `test_methods.R`.

## State
- Uncommitted: clean. Local-only branches: `feat/model-benchmark` @ `87181b5`, `feat/sofa-comparator` @ `7efc410` (no PR).
- Tests: 155 expectations, 0 failures: `Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'` with `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin/Rscript` (no `Seurat` in that env). Toy methods check: `R_LIBS_USER=/nonexistent Rscript benchmarks/coordination/test_methods.R`.
- Running: nothing here. Biomni session `sess_a9d08850c561` (project `prj_0135m4eKKYB2YYYOelt1y9Xo`) is paused at its 750-credit limit ("send a message to continue"); my Part 2-4 request is unanswered. Its agent made branch `feat/coordinated-programs` (57 files on `b0c6986`, push got 403); not applied here, review before use.
- `.c1-work/` (295 MB) is a separate clone from another session with unpushed C1 commits `0ffa453`, `e95f1d7` and uncommitted edits. Hidden via `.git/info/exclude`. **Do not delete.** Same session saved `reports/plans/2026-10-09-selective-coupled-ebmf-gene-sharing.md` (not reviewed).
- SLURM `--partition=all`; MaxArraySize 125, so arrays go in chunks (`OFFSET` in `run.sbatch`, D7).

## Locked decisions
- `strata` rejected as registered, no post-hoc guard (user, 2026-10-10; D8). Future guards relative to `base`.
- Metaprogram path: repair in-repo methods first, then register; `ebmf2` and external comparators later (user, 2026-10-10).
- Biomni files: user downloads them; no more Biomni credits spent for now (user, 2026-10-11).
- Coverage guardrail `frac_untested <= 0.05`, F as registered (user, 2026-10-08).
- Main run dropped `k_stab` and `w_S` (runtime guard, D6) (user, 2026-10-08).
- Labels `coordinated`, `gene_reuse`, `both`, reported separately (user, 2026-10-09).
- HVGs Seurat v3 `vst`, union of per-cell-type HVGs; test route A (fit on union) vs B (refit loadings on fixed scores) before choosing (user, 2026-10-09).
- No MOFA-on-scores or GFA; MOFAcellulaR only as optional comparator (user, 2026-10-09).
- Defaults: `covariate_mode="residualize"`, K=10, `point_laplace`, `var_type=1`; keep `fit_celltype_programs`; L1/L2 priors and EV-BIDIFAC rejected.

## Findings to carry (toy, 10 seeds, 999 perms, recovery R100/R50/R20)
- `node_cor` 1/1/1, `gene_match` 1/1/1, `pa_cc` 1/0/0 (merges everything), `subspace_sum` after deflation 1/1/0.3. `subspace_sum` is the weak arm; cheaper to demote than tune.
- Biomni finalist `n_agg_noanc` (rank-Gauss then residualize on strata; ambient filter; node statistic `A_v = sum_c' (max_u |T_vu| - 2)_+^2`; max over nodes under within-strata permutation) targets exactly R50/R20 power. Mixing edge and node test families gave null FWER 0.14. Their scoring rule (|cor| >= 0.4, SITE pseudo-signal, F1 = max over CMPs): `reports/biomni_cmp_lessons.md` B1b.
- Their real data at cap 30: new CMPs are mostly donor-genotype (GSTM1, HLA, CCL3L3, sex chr) or Ig/TCR, not coordinated programs; their simulator misses both (B1c). Stated by them, not re-run.

## Dead ends — do not redo
- `p_pair` alone to compare runs (saturates); AUROC of `z_pair` (1.00 unadjusted); use `excess_frac`/`z_pair`, FPR, power.
- `ebnm_flat` for fixed covariate factors (NA ELBO); use fixed-scale `ebnm_normal`.
- A floor on `arm_k_stab`.
- Complete-case methods (DIALOGUE-style multi-CCA, ICA, RaJIVE) as main route: only 9 of 124 COMBAT donors have all 10 cell types (`biomni/REVIEW.md`).
- CP/Tucker tensor as main route; CRAN `multiCCA`, `mbpls`.
- Gaussian-z simulator with ICA: unidentifiable; plant non-Gaussian z if ICA is benchmarked.
- Messaging Biomni while it is paused (402, then credit-limit pause).

## Read first
- `reports/biomni_cmp_lessons.md` — Part A ours, Part B theirs (B1b callers and scoring, B1c real data), adoption list in Part C.
- `reports/biomni_transfer/part1_callers_scoring.md` — their Part 1 reply verbatim, with R code.
- `reports/metaprogram_methods_review.md` — §5 phase 1 plan, §6 open decisions, §10 HVG test.
- `benchmarks/coordination/PREREG.md` — draft, four blockers: method list vs `methods.R`, replicates and paired inference, node-to-truth rule (their rule can fill it), pinned EBMF config.
- `benchmarks/coordination/methods.R:156` (`subspace_sum`), `:81` (`node_cor`), `:26` (`.cluster_edges`), `:59` (`.maxT_edges`).
- `benchmarks/model_selection/DEVIATIONS.md` D1-D8, `summarize.R:65-76` (untested guard).

## Open items
- Phase 1 remaining: cluster-level coherence test (candidate: `node_agg`), node density filter, PREREG blockers, then register `prereg-coord-v1`.
- Adoption items not needing Biomni data: untested guard relative to `base` in PREREG; Ig/TCR-like and donor-deletion-like signals in `sim_coord.R`; report K-cap hits in real runs; genotype-marker flag in the real-data CMP table; never change two factors between compared runs.
- Toy finding to verify: stage-2 EBMF (`ebmf2`) on stacked node scores found a 3-cell-type program (10 of 12 seeds) but not a 2-cell-type one with 22% missing donors.
- Count-level simulator (current is Gaussian); install Seurat; check whether `benchmarks/biomni_replication/data/combat/pb_counts.parquet` is raw counts.
- `both` label merge rule; scoring tiers for methods without node output; strata null for external arms.
- Biomni M4 pre-work if a fit change is ported: `flash_control` conflict detection (`R/ebmf.R:123-132`), package version in `provenance`, small-n K guard.
- Biomni M2 findings not acted on: confounder check not in `summarize.R`; NA passes some guards; no multiplicity control across arms.
- `feat/sofa-comparator`: open a PR or archive it.
