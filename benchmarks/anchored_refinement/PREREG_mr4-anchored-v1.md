# Pre-registration: MR-4 anchored refinement + ablation (REG-MR4-anchored-v1)

Registered in `reports/STAGE_B_ROADMAP.md` (REG-MR4-anchored-v1); this file
restates the rule for the run harness and freezes the design. Any deviation
goes in `DEVIATIONS.md` next to it. Written 2026-10-22 before any run.

## Question

Does anchored refinement — initialize Stage A programs and refit on
cell-level (or metacell-level) counts with `flashier.utils::fit_initialized_ebmf()`
— improve program recovery over the de-novo cell-level fit, without drifting
the programs or degrading reconstruction?

## Anchor

Stage A v0.2.0 operating point: donor pseudobulk per cell type,
`fit_celltype_programs()` with `max_factors = 30` (the cmp-bench-v1.2 winner
cap), `loading_prior = "point_laplace"`, `var_type = 1`, 2000 variable genes,
`covariates = "site"`, `covariate_mode = "residualize"`, `seed = 1`, then
`canonicalize_programs()`. Cell-level activities are initialized with
`project_programs_to_cells()` (the package's projection oracle).

## Arms (all paired on the identical cell-level counts)

| Arm | Fit |
|---|---|
| `denovo_cell` | `fit_celltype_programs()` on cell-level data (comparator; the Stage B de-novo operating point) |
| `anch_free` | `fit_initialized_ebmf(fixed_side = "none", add_factors = 0)` |
| `anch_fixed_act` | `fixed_side = "activity"` (anchor activities fixed, loadings refit) |
| `anch_fixed_eff` | `fixed_side = "effects"` (anchor loadings fixed, activities refit) |
| `anch_add5` | `fixed_side = "none", add_factors = 5` |
| `anch_free_mc` | metacell-level counts (~10 cells/metacell, k-means grouting; SuperCell on the real run), `fixed_side = "none", add_factors = 0`, `S = sqrt(size)` |

All anchored arms: `effect_prior = ebnm::ebnm_point_laplace`,
`activity_prior = ebnm::ebnm_point_normal`, `var_type = 1`,
`feature_scale = "logcounts"`, `seed = 1`.

## Data

Semi-synthetic cell-level counts (the harness simulator): 4 cell types
(A–D), 60 donors (2 sites), G genes; a 50-gene shared program planted in the
A–B pair with donor activities z ~ N(0,1); 2 private programs per type;
cells per donor lognormal (median 200, sd_log 0.8); Poisson counts on
exp(log-profile). Amplitude chosen by pilot so `denovo_cell` CMP power is
~0.5 (stored in `amplitudes.csv` before the main run; pilot reps discarded).
The real-data variant (anchor = Stage A fits on real pseudobulk, query = the
cohort's single-cell counts) uses the same harness with the simulator
swapped; CMP F1 is then replaced by the drift + stability + reconstruction
endpoints (no planted truth on real data).

## Replicates

200 replicates, seeds `4400000 + rep` (pilot 4390001–4390020, discarded).
Every arm runs on the identical dataset in each replicate (paired).

## Metrics

- **Held-out CMP F1** (primary): `find_coordinated_programs()` (the
  `n_agg_noanc` caller settings: rank-first, ambient "none" for parity across
  arms, `n_perm = 199`, `alpha = 0.05`, `min_shared = 20`, `min_donors = 20`;
  lowered proportionally on the pilot-scale smoke) on each arm's
  scores/loadings wrapped in the fit contract; F1 over the 6 unordered pairs
  against the planted A–B pair.
- **Drift**: per cell type, Hungarian-matched mean |cor(anchor donor scores,
  refit donor-mean activities)| over anchor programs; reported as
  `1 - value`.
- **Per-cell activity stability**: 10 subsamples of 80% of cells per type;
  re-project the arm's loadings onto the subsample; mean Spearman of
  per-cell activities on shared cells over programs and subsamples.
- **Reconstruction** (guardrail): per cell type, mean per-gene residual
  variance ratio vs `denovo_cell` on the shared gene panel (ELBO-equivalent
  proxy; degradation > 1% flags the guardrail).

## Decision rule

Adopt anchored refinement into `R/` iff, vs `denovo_cell`:

1. held-out CMP F1 difference >= **+0.02 absolute**, AND
2. paired bootstrap 95% CI (2000 resamples) of the F1 difference excludes 0, AND
3. no reconstruction degradation > 1% (residual-variance ratio <= 1.01) in
   either cell type, AND
4. drift is reported but not gated (measurement endpoint).

Ablation arm (factor-utility scoring from the model_selection factorial —
`k_stab`-style stability filter at `recovery_freq >= 0.7`, rerun on the
Stage B operating point) is evaluated with the same rule and reported
separately; it does not gate the anchored-refinement decision.

## Cost

Dominated by the cell-level refits; estimate from the model_selection main
(~324 CPU-h) scaled by cell count — budget ~400–600 CPU-h for both arms on
the Slurm cluster. **Staged scripts, not launched from the sandbox**; the
in-sandbox deliverable is the pilot-scale smoke.

## Operationalization note (2026-10-22, pre-run; smoke-driven, no rule change)

recon_ratio is computed on the gene panel SHARED between the arm's loadings
and the de-novo fit's variable panel (the two panels are selected on different
variance bases: anchor pseudobulk vs cell-level). Panels are never identical,
and character row-indexing errors on missing names, so the shared-panel
intersection is the only well-defined version of the registered
"per-gene residual variance ratio vs de-novo" guardrail. If the shared panel
is < 10 genes the ratio is NA (excluded from the min). Drift's Hungarian
matching is over anchor programs (K = min(dim) after column filtering).
