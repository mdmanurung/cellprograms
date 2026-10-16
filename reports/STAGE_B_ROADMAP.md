# Stage B Roadmap: from donor pseudobulk programs to cell-level programs

Status as of 2026-10-16. Stage A (donor-pseudobulk EBMF backend, canonicalization,
stability, CMP detection, coordination arms) is feature-complete and released through
v0.2.0. Stage B extends the fitted programs to the cell level and hardens the
methodological story. Each milestone below follows the repo law: the registered rule
is written before any R/ change, and only rule-clearing results enter the package.

## MR-0 — Frozen preprocessing specifications (IMPLEMENTED this phase)

`freeze_preprocessing()` resolves the feature panel and adjustment decisions once and
returns a content-addressed spec (`cell_program_preproc_spec`, sha256 digest).
`fit_celltype_programs(preprocessing_spec = ...)` replays it fail-closed: conflicting
explicit arguments are refused, spec cell types must cover the dataset, and all frozen
genes must be present. Refits are bit-identical to explicit-argument fits (tested).
The spec digest is recorded in `fit$preprocessing$frozen_spec`.

## MR-1 — Pseudobulk-to-cell projection (IMPLEMENTED this phase)

`project_programs_to_cells()` puts single-cell query data into a fit's feature space
by replaying the recorded preprocessing (per-gene OLS residualization coefficients,
then stored centering/scale), infers activities against the frozen loadings with
`flashier.utils::project_factor_activity()`, and canonicalizes with the same RMS/sign
multipliers as the donor scores. Design decisions (all fail closed):

- Requires a canonicalized fit with a flashier.utils manifest and declared
  `feature_scale`; legacy-backend and fixed-covariate (`covariate_mode = "fixed"`)
  fits are refused — their adjustment cannot be replayed.
- Covariate replay rebuilds the design from `donors` + `sample_metadata` and requires
  the design columns to match the fit's stored columns exactly (a level lost or gained
  among query donors aborts rather than silently projecting through a different
  adjustment).
- Query genes are matched by intersection with `min_coverage` (default 0.9); the
  oracle's rank/coverage guards apply unchanged.
- Output: cells x programs activity matrix (NA off-cell-type), per-type
  reconstruction metrics, coverage table, analysis metadata, provenance.

Verified: projecting the training pseudobulk reproduces canonical donor scores
(cor > 0.99 per program); wrapper output is bit-identical to manual replay +
direct oracle + manual canonicalization.

## MR-2 — cNMF baseline (STAGED)

Registered comparison: cNMF (cottrill & smibert) on the same pseudobulk inputs,
matched K grid, same donor scores -> cell projection path as the flashier backend.
Rule: cNMF replaces the default backend only if it beats the registered Stage A
operating point on the cmp-bench-v1.2 held-out CMP F1 by >= 2 points absolute with
non-overlapping bootstrap CIs. Cost: cNMF is cheap (~minutes per dataset on 8 cores);
benchmark harness ~15 CPU-h. Registration draft to be filed before any benchmark run.

## MR-3 — Metacells design doc (STAGED)

Design-only milestone: when single-cell factorization is needed (sparse cell types,
rare states), metacells are the candidate middle ground between donor pseudobulk
(loses within-donor heterogeneity) and raw cells (noise). Deliverable: a design doc
covering metacell grouting (SEACell/Metacell-2 candidates), program fitting on
metacell pseudobulk, and projection of remaining cells; no code until the doc is
reviewed. Cost: authoring only.

## MR-4 — mcRigor / anchored refinement + ablation benchmarking (STAGED)

Two registered arms, both built on flashier.utils `fit_initialized_ebmf()`:
(a) anchored refinement — initialize Stage A programs, refit on cell-level or
metacell-level data with programs fixed/anchored, measure drift; (b) ablation
benchmarking — the OPT utility consolidation (factor-utility scoring arms from the
model_selection factorial) rerun on the Stage B operating point. Rule for adopting
anchored refinement into R/: held-out CMP F1 improvement >= 2 points absolute,
bootstrap CI excluding 0, and no reconstruction degradation > 1% ELBO-equivalent.
Cost: dominated by the cell-level refits; estimate from the model_selection main
(~324 CPU-h) scaled by cell count — budget ~400-600 CPU-h for both arms on the
Slurm cluster (staged scripts, not launched from the sandbox).

## Cost summary for the staged main runs

| Run | Est. cost | Status |
|---|---|---|
| prereg-coord-v1 main (coordination arms) | ~5-6 h wall, ~85 s/fit at G=1000 pilot scale | staged (pilot.R ready) |
| MR-2 cNMF benchmark | ~15 CPU-h | registration draft below |
| MR-4 anchored + ablation | ~400-600 CPU-h | staged |

## Registration drafts

- **REG-MR2-cnmf-baseline-v1**: datasets = COMBAT pseudobulk (sex contrast) +
  Stephenson pseudobulk; K grid = {5, 10, 20, 30}; metrics = CMP membership F1
  (held-out), reconstruction, stability (existing arm metrics); decision rule above.
- **REG-MR4-anchored-v1**: anchor source = Stage A v0.2.0 fits; query = single-cell
  (or metacell) counts; metrics = drift (canonical score correlation), held-out CMP
  F1, per-cell activity stability under cell subsampling; decision rule above.

## Gates

- Stage B code milestones are gated behind the Stage A release (v0.2.0 tag).
- Every benchmark-gated change needs its registered rule in this file (or a successor
  commit) BEFORE the R/ change lands.
- Nothing red merges: full testthat suite + R CMD check must pass on every push.
