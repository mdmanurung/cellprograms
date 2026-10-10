# Handoff — cellprograms: strata rejected; metaprogram work = repair in-repo methods first

**Date:** 2026-10-10 · **Branch:** `feat/coordination` (`origin/main` is at `b0c6986`, 3 commits behind; fix is `10c9b55`) · **Status:** M0, M1, fix plan, M2 done, no arm adopted. Both open decisions answered (below). Repair phase started.

## Goal
(1) Decide which fit modifications earn their complexity (model-selection benchmark, pre-registered at `prereg-v1`). (2) Land a defensible way to find metaprograms (programs shared, coordinated or recurring in all or only some cell types) from separate per-cell-type EBMFs. Nothing enters `R/` unless it clears a registered rule.

## Next action
Work `docs/metaprogram_methods_review.md` §5 phase 1 in order. Done: `subspace_sum` deflated sequential test (`10c9b55`). Next: (a) decide what to do about `subspace_sum` R20 power (3/10 at 999 perms, see Open items); (b) cluster-level coherence statistic + test; (c) cNMF-style node density filter; (d) the four `PREREG.md` blockers (method list vs `methods.R`, replicates and paired inference, node-to-truth matching, pinned EBMF config); then register `prereg-coord-v1`.

## State
- Uncommitted: clean after this handoff commit. Local-only branches: `feat/model-benchmark` @ `87181b5`, `feat/sofa-comparator` @ `7efc410` (no PR).
- Tests: 155 expectations, 0 failures: `Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'` with `/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/cellprograms-r/bin/Rscript` (no `Seurat` in that env).
- Running: nothing. SLURM `--partition=all`; MaxArraySize 125, so arrays go in chunks (`OFFSET` in `run.sbatch`, D7).
- `.c1-work/` (295 MB) is a separate clone from another session with unpushed C1 commits `0ffa453`, `e95f1d7` and uncommitted edits to `R/coordination.R`, `README.md`, `benchmarks/coordination/PREREG.md`. Hidden locally via `.git/info/exclude`. **Do not delete.**
- That session also saved `docs/plans/2026-10-09-selective-coupled-ebmf-gene-sharing.md` (gene-loading recurrence workstream; next step GS-00 spec only, not preregistered). Not reviewed by me.

## M2 results (`benchmarks/model_selection/RESULTS.md`, `results/decisions*.csv`)
- 900 datasets, 8 arms, no failures. As registered: no arm adopted, M3 not triggered.
- `k_n` dP on S3 +0.20 `param` / +0.50 `semi`, but `semi` S0 F 0.14 -> 0.705; `k_20`, `genes_after` no effect; `vt_12` harmful (dR -0.41 / -0.54).
- `strata` dF on S0 -0.24 `param` / -0.075 `semi`; blocked only by the `frac_untested <= 0.05` guard, which fails every `param` arm (12-27% untested, `base` included).
- `base` F on S0: 0.295 `param`, 0.14 `semi`. Confounder check: `none` F = 1.0 on `param` (passes), 0.09 on `semi` (fails; weak vs unaligned institute effect not distinguishable).

## Locked decisions
- `strata` rejected as registered, no post-hoc guard (user, 2026-10-10; D8 in `DEVIATIONS.md`). Future benchmarks should register the untested guard relative to `base`.
- Metaprogram path: repair in-repo methods first, then register; stage-2 EBMF (`ebmf2`) and external comparators come after (user, 2026-10-10).
- Coverage guardrail `frac_untested <= 0.05`, F left as registered (user, 2026-10-08). Known flaw: absolute, not relative to `base`.
- Main run dropped `k_stab` (9.7x time) and `w_S` (4.5x) (user, 2026-10-08).
- Metaprogram labels `coordinated`, `gene_reuse`, `both`, reported separately; gene reuse is an accepted target (user, 2026-10-09).
- HVGs: Seurat v3 `vst`, union of per-cell-type HVGs; test route A (fit on union) vs B (refit loadings on fixed scores) empirically before choosing (user, 2026-10-09).
- No MOFA-on-scores or GFA; direct MOFAcellulaR only as an optional power/null comparator (user, 2026-10-09).
- Hard adjustment (`covariate_mode="residualize"`) default; repo defaults K=10, `point_laplace`, `var_type=1`; do not replace `fit_celltype_programs`; L1/L2 priors and EV-BIDIFAC stay rejected (ICA/RaJIVE/jNMF allowed as comparator arms only).

## Dead ends — do not redo
- `p_pair` alone to compare runs (saturates); AUROC of `z_pair` (1.00 even unadjusted); use `excess_frac`/`z_pair`, FPR and power.
- `ebnm_flat` for fixed covariate factors (NA ELBO); use fixed-scale `ebnm_normal`.
- A floor on `arm_k_stab` to rescue it.
- Complete-case methods (DIALOGUE-style multi-CCA, ICA, RaJIVE) as the main route: only 9 of 124 COMBAT donors have all 10 cell types (`biomni/REVIEW.md`).
- CP/Tucker tensor as the main route: assumes the same genes and activity in every member cell type; CRAN `multiCCA` (kernel/functional) and `mbpls` (supervised) do not fit.
- Gaussian-z simulator with ICA: unidentifiable; plant non-Gaussian z if ICA is benchmarked.

## Read first
- `docs/biomni_cmp_lessons.md` — Biomni CMP benchmark lessons; Part A (ours) done, Part B (theirs) waits on a Biomni retry scheduled 2026-10-11 00:30 (session-only cron `7e52e47d`; first send got HTTP 402) or on downloaded `tusoskill/cmp_v1`, `cmp_v2` files.
- `docs/metaprogram_methods_review.md` — candidates, per-framework verdicts, staged plan and stopping rule (§5), open design decisions (§6), upstream HVG test (§10).
- `benchmarks/coordination/PREREG.md` — draft, unregistered; four blockers: method list vs `methods.R`, replicates and paired inference, node-to-truth matching rule, pinned EBMF config.
- `benchmarks/coordination/methods.R:146-195` (`subspace_sum`; `cummax` + `break` is not a valid step-down), `methods.R:11-17` (`.perm_rows` fix).
- `benchmarks/model_selection/DEVIATIONS.md` D3-D7, `summarize.R:65-76` (untested guard).

## Open items
- `subspace_sum` after deflation (10 toy seeds, 999 perms; script was in scratchpad, not kept): R100 and R50 recovered 10/10, no wrong members, null gives 0 clusters; R20 (2 cell types) recovered only 3/10. The earlier `ct01,ct09` was a 199-perm one-off. Same toy data, 999 perms, recovery R100/R50/R20: `node_cor` 1/1/1, `gene_match` 1/1/1, `pa_cc` 1/0/0 (merges everything), `subspace_sum` 1/1/0.3. So R20 is detectable and `subspace_sum` is the weak arm (max-eigenvalue null), not the data. Cheaper to drop or demote it than to tune it.
- Toy finding to verify: stage-2 EBMF on stacked node scores found a 3-cell-type program (10 of 12 seeds) but not a 2-cell-type one with 22% missing donors; complete donors found both (about 90%). Needs a real benchmark arm (`ebmf2`).
- Count-level simulator (current one is Gaussian), install Seurat, check whether `benchmarks/biomni_replication/data/combat/pb_counts.parquet` is raw counts; HVG before vs after adjustment depends on the `genes_after` result (no effect, so before).
- `both` label merge rule (groups intersect / contain / overlap threshold); scoring tiers for methods without node output; strata null for external arms.
- Biomni M4 pre-work if any fit change is ported: `flash_control` vs argument conflict detection (`R/ebmf.R:123-132`), package version in `provenance`, small-n behaviour of a K <= n-q-1 guard (error vs cap).
- Biomni M2 findings not acted on: confounder check not computed by `summarize.R`; NA passes some guards; no multiplicity control across arms.
- `feat/sofa-comparator`: open a PR or archive it.
- From the Biomni CMP benchmark session `sess_a9d08850c561` (sibling effort; I read the final reply and file names only, not the files): (1) K cap 10 truncated real data (4 cell types hit 40 at cap 40; CMPs 13 -> 31) while the simulator screen was flat from 30 to 40, so a sim-flat `k_*` result does not show real data is untruncated; (2) IGHV/IGKV/IGLV/TRBV-dominated programs rose 4 -> 6 (Stephenson) and 1 -> 9 (COMBAT) at the higher cap and the simulator has no such signal: filter them before HVG selection as a separate labelled arm, and check the Seurat `vst` union does not feed them into the gene-reuse label; (3) their dev / held-out seed split and a real-data CMP table (ambient, Ig/TCR, severity, site replication) are worth copying for our real-data stage; (4) never compare runs that changed two things (their COMBAT cap-10 vs cap-40 also changed sex adjustment).
