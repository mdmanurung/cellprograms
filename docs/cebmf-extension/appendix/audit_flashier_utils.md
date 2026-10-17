# Audit: flashier / flashier-utils / cellprograms — seams for an alternative EBMF backend

Sources inspected (read-only):
- flashier source: `/workspace/flashier_src_audit` (GitHub stephenslab/flashier, master @ 3969ccd, DESCRIPTION Version 1.0.59; installed lib is 1.0.60 at `/workspace/.Rlib/flashier` — lazy-loaded, no R sources on disk; patch-level difference only, public API identical). **Note: `flash_r1()` / `flash_r2()` do NOT exist in flashier 1.0.x** — they were removed after 0.6-xx; grep of `NAMESPACE`/`R/` finds no such export. The 1.0 API is `flash()` + pipeable `flash_xxx()` verbs.
- flashier-utils: `/workspace/flashier-utils_audit` (cloned, HEAD = 271f583 "docs(reports): GO memo…", Version 1.1.0.9000).
- cellprograms: `/workspace/repo2_full` (Version 0.2.0, untouched).

---

## (a) flashier API + fit-object inventory

### Public API
| Item | Location |
|---|---|
| `flash(data, S=NULL, ebnm_fn=ebnm_point_normal, var_type=0L, greedy_Kmax=50L, backfit=FALSE, nullcheck=TRUE, verbose=1L)` | `R/flash.R:216-241` (docs 1–214) |
| `flash_init(data, S=NULL, var_type=0L, S_dim=NULL)` | `R/flash_init.R:33-63` |
| `flash_r1` / `flash_r2` | **removed in 1.0** (not in `NAMESPACE:35-80`); replaced by `flash()` + `flash_greedy`/`flash_backfit` pipeline |
| Greedy/backfit/nullcheck verbs | `R/flash_greedy.R`, `R/backfit.R:1-88`, `R/flash_nullcheck.R:36-105` |
| `flash_ebnm(...)` helper (wraps `ebnm::ebnm`/`ebnm_group`; custom fns allowed) | `R/flash_ebnm.R:44-100` (`handle.ebnm.fn` 102–117, probe `test.ebnm.fn` 121+) |
| Accessors: `flash_fit_get_pm/p2m/est_tau/fixed_tau/tau/elbo/KL/g` | `R/flash_fit_object.R:47-88` |

### `flash` object fields (built in `wrapup.flash`, `R/wrapup.R:1-95`)
- `n_factors` (L16), `pve` (L15, `calc.pve` L111+), `elbo` (L18) — **single final value, not a trace**; the per-iteration trace lives only in verbose output (`R/flash_verbose_fns.R`).
- `residuals_sd` (L20-25; scalar `1/sqrt(tau)` if `is.tau.simple`, else per-mode list).
- `L_pm, F_pm` (L32-33), `L_psd, F_psd` (L34-35, `sqrt(pmax(EF2-EF^2,0))`), `L_lfsr, F_lfsr` (L36-44, only if `all.valid`), `L_ghat, F_ghat` fitted priors (L67-71), `sampler` (L78), `flash_fit` (L88, internal `flash_fit` object with `EF`, `EF2`, `KL`, `g`, `tau`, `est.tau`, `given.tau`, `R`, `Y`, `Z`/nonmissing, `fix.dim`/`fix.idx`).
- "fixed elements": `flash$fix.dim`/`flash$fix.idx` set in `flash_init.R:52-53`; fixed factors via `flash_factors_fix` (used by cellprograms `.flash_fixed`, `R/ebmf.R:429-455`).

### var_type → noise structure
Validated in `R/arg_checks.R:61-69` (`must.be.valid.var.type`: integers 0..dim, no repeats, 0 must be alone). Semantics (doc `R/flash.R:139-150`):
- `NULL`: `S` accounts for *all* residual variance (nothing estimated).
- `0`: single constant variance added to `S` (default).
- `1`: one variance per row (donor); `2`: one per column (gene); `c(1,2)`: Kronecker rank-one `s_ij = a_i*b_j` (slow).
- Storage decision `is.S2.stored` in `R/set_data.R:118-124`; Kronecker + missing data unsupported (`R/set_data.R:40-44`).

### Missing data
Handled in `set.flash.data` (`R/set_data.R:20-25`): `Y[is.na(data)] <- 0`, `Z <- 1*!is.na(data)`; weighted updates via `Z` throughout. `flash_init` itself has no NA branch — it delegates to `set.flash.data` and precomputes `n.nonmissing` (`flash_init.R:66-80`). All-zero rows/cols rejected unless `S>0` or `var_type=0` (`arg_checks.R:71-94`). Low-rank `u/d/v` inputs must be NA-free (`set_data.R:103-105`).

### Where `S` (per-gene variance floor) enters
`S` is an argument of `flash()` (`flash.R:216`) and `flash_init()` (`flash_init.R:33`), consumed only by `set.flash.data` (`set_data.R:7-60`): scalar/vector/matrix → `given.S2` (when `is.S2.stored`) or `given.tau = 1/S^2` (`set_data.R:45-58`); vector `S` aligned via `infer.S.dim` (`set_data.R:96-110`), square matrices need `S_dim`. **cellprograms never passes `S`** — `args` in `R/ebmf.R:289-297` contains no `S`, so all residual variance is estimated by flashier (`var_type=1` default). A per-gene variance floor would enter a Python backend as a per-gene `tau` floor, i.e. the `given.tau`/`S` channel.

### Backfit / strong rule / nullcheck
- Backfit: cyclic `update.one.factor` per k (`R/backfit.R:12-29`), accept only if objective improves/valid (`backfit.R:17-24`); `update.factors.in.kset` recomputes R, tau, obj (`backfit.R:1-11`).
- Nullcheck ("strong rule"): `flash_nullcheck` (`R/flash_nullcheck.R:36-81`) zeroes each factor (`nullcheck.factor` L85-99) and removes it if ELBO does not drop by more than `tol = sqrt(.Machine$double.eps)*prod(dims)`; `flash_factors_remove` re-runs `nullcheck.factor(..., tol=Inf)` (`R/flash_factors_remove.R:68`).

### ebnm prior families
Any `ebnm::ebnm_*` function works as `ebnm_fn` (point_normal default, `flash.R:216`; point_laplace, point_exponential, unimodal, unimodal_nonnegative, normal_scale_mixture, horseshoe, semiconjugate, flat, etc. — everything exported by package `ebnm`). Per-side and per-k specification via a length-2 list (`handle.ebnm.fn`, `flash_ebnm.R:102-117`). flashier-utils records the matched `ebnm` export name or `"custom_callback_required"` (`fit.R:139-144`).

---

## (b) flashier-utils exports and role (HEAD 271f583)

flashier-utils is a **post-fit analysis + provenance layer over flashier** (Imports: flashier, ebnm, Matrix, digest, ggplot2). It does not reimplement EBMF; it wraps `flashier::flash*` calls, attaches a `flashbridge_metadata` attribute to the returned flash object, and adds single-cell-oriented tooling. Key groups (file: def line):

**Fitting wrappers** (all return a native flash object + metadata attr):
- `fit_ebmf` (`R/fit.R:80-171`): validates matrix/ids, forces explicit `sample_side`, optional feature-side intercept (`flash_add_intercept`, fit.R:113-131), seeds via `.with_seed`, then either `flashier::flash(...)` (fit.R:126-128) or the decomposed `flash_init→flash_greedy→flash_backfit→flash_nullcheck` pipeline with `control$greedy/backfit/nullcheck` (fit.R:129-140, 108-111). Attaches `flashbridge_metadata`: `fit_id` (sha256), `factor_ids` ("F1..Fk"), `sample_side`, row/column ids, `feature_scale`, `prior_support` (activity/effects), `prior_identifiers`, seed + `RNGkind()`, full `configuration` (fit.R:135-158).
- `fit_multistart_ebmf` (`R/fit.R:228+`): seeded candidate starts (`native_random`/`softImpute`/`irlba` initializers, `.multistart_initializer` fit.R:173-197), ELBO-ranked selection.
- `fit_initialized_ebmf` (`R/initialized.R:37`), `fit_nonnegative_ebmf` / `fit_seminonnegative_ebmf` (`R/operations.R:99,134`), `refit_factors`, `fix_factors`, `remove_factors`, `reorder_factors` (`R/operations.R:182,219,264,295`).

**Extraction / normalization:** `factor_activity`, `factor_effects`, `factor_uncertainty`, `factor_pve` (`R/extract.R:52,78,112,169`); `rank_features`, `factor_features`, `factor_distinctiveness` (`R/features.R:30,130,73`); `standardize_factors` (`R/representation.R:183`; scalings `ldf`/`rms`, orientation `anchor_feature` — used by cellprograms canonicalization).

**Projection (single-cell):** `project_factor_activity` (`R/projection.R:229+`) and `project_factor_features` (`R/projection.R:152`): regress query genes on the reference effect matrix (`B = source$effects[,k]`), rank-deficiency/coverage gates (`qr(shared)$rank<ncol(shared)` → stop), memory budget, optional `S` and `var_type` for the query fit; returns `representation$activity`, `reconstruction` metrics, `analysis_metadata`. `backproject_contrast` (`R/backprojection.R:48`), `validate_backprojection`.

**Diagnostics/annotation/stats:** `factor_holdout` (`R/holdout.R:86`, native projection CV), `factor_stability`, `factor_replicability` (`R/alignment.R:115`), `factor_cor` (`R/correlation.R:32`), `test_factors` (limma/varPartition designs, `R/association.R:101`), `summarize_factor_activity`, `rank_factors_by_metadata`, `compare_factor_groups` (`R/metadata.R:64,146,200`), `enrich_factors` (`R/enrichment.R:78`), `factor_network`, `factor_diagnostics`, `factor_specificity`, `factor_intervals`, `factor_gate_alignment`, `embed_factors`, `shrink_factor_effects` (`R/shrinkage.R:38`), `provenance` (`R/provenance.R:19`), ~15 `plot_*` functions (`R/plots.R`, `R/programs.R`, `R/program-markers.R`, `R/program-diagnostics.R`). Full export list: `NAMESPACE:3-55` (55 exports).

---

## (c) cellprograms seam map for backend substitution

Pipeline (`/workspace/repo2_full`):

1. **Input/normalization** — `fit_celltype_programs` (`R/ebmf.R:94-369`): per cell type, subset genes (`.select_features` `ebmf.R:497-527`), optional OLS residualization (`ebmf.R:247-253`, `.cov_design` `ebmf.R:381-427`), center/scale with replay record (`ebmf.R:255-263`: `replay = list(center, scale, residualize_coefficients, design_columns)`), drop zero columns (`ebmf.R:266-281`), C4 residual-df guard (`ebmf.R:284-299`).
2. **Backend dispatch — THE SEAM** (`ebmf.R:289-293`):
   ```r
   ordinary <- is.null(Xc) && use_fu          # ebmf.R:289
   run_flash <- function(args) {              # ebmf.R:290
     if (!is.null(Xc)) return(.flash_fixed(args, Xc))
     if (ordinary) return(.fit_via_flashier_utils(args, seed, feature_scale))
     do.call(flashier::flash, args)           # ebmf.R:293 legacy backend
   }
   ```
   `use_fu` gate: `ebmf.R:126-127` (`getOption("cellprograms.use_flashier_utils", TRUE) && requireNamespace("flashier.utils")`). A Python backend replaces the body of `run_flash` (or adds a third branch before the flashier ones).
3. **flashier call sites:** legacy `do.call(flashier::flash, args)` `ebmf.R:293` with `args` built at `ebmf.R:289-297` (`data=Y, ebnm_fn=list(ebnm_normal, prior_fn), var_type, greedy_Kmax, backfit, nullcheck, verbose=0`); flashier.utils route `.fit_via_flashier_utils` `ebmf.R:467-495` → `flashier.utils::fit_ebmf(X=args$data, sample_side="rows", ebnm_fn, var_type, max_factors, backfit, nullcheck, seed, verbose, feature_scale)` at `ebmf.R:479-489`; fixed-covariate path `.flash_fixed` `ebmf.R:429-455` (`flash_init | flash_set_verbose | flash_factors_init | flash_factors_fix | flash_greedy` at `ebmf.R:432-438`, then backfit/nullcheck `ebmf.R:440-442`).
4. **Fallback ladder** (`ebmf.R:300-338`): attempts = [prior as requested] → [point_normal if point_laplace] → [each ± no_nullcheck]; each attempt re-seeded (`set.seed(seed)` `ebmf.R:313,323`), failure = error **or** non-finite ELBO via `.fit_ok` (`ebmf.R:374`); winner recorded in `fits[[ct]]$fallback` (`ebmf.R:336-341`). A Python backend must either honor this ladder or map its own failure modes onto `fallback` labels.
5. **Canonicalization** — `canonicalize_programs` (`ebmf.R:545-627`): reads `fit$fits[[ct]]$flash$L_pm/F_pm` (`ebmf.R:556-558`); if `attr(f,"flashbridge_metadata")` present routes through `flashier.utils::standardize_factors(f, sample_side="rows", scaling="rms", orientation="anchor_feature")` (`ebmf.R:561-570`), else inline RMS `s_k = sqrt(colMeans(Z^2))` + anchor-gene sign (`ebmf.R:571-587`); writes `fit$canonicalization[[ct]] = list(Z, W, scale_multiplier, sign_multiplier)` (`ebmf.R:596-600`) and `fit$scores[[ct]]` (obs×K, rownames=observation ids, colnames `paste0(ct,"_",1:K)`) / `fit$loadings[[ct]]` (genes×K) (`ebmf.R:602-625`).
6. **Projection/scoring** — `project_programs_to_cells` (`R/cell-projection.R:80-227`): **hard-requires** `fit$scores` + `fit$canonicalization` (`cell-projection.R:88-90`) and `attr(fit$fits[[ct]]$flash, "flashbridge_metadata")` + `preproc_replay` (`cell-projection.R:141-152`); replays residualize→center→scale on the query (`cell-projection.R:160-175`), then calls `flashier.utils::project_factor_activity` (`cell-projection.R:177-186`) and re-applies canonicalization multipliers (`cell-projection.R:188-191`). Other consumers: `program_scores/program_loadings/program_summary/top_genes` (`R/extract.R:37,67,88,118`), CMP interpretation `annotate_cmp_tf`/`.cmp_oriented_loading` (`R/cmp-interpret.R:144,49-54`), annotation `annotate_cmp_lr/themes`, `cmp_triage` (`R/cmp-annotate.R:80,208,321`), coordination `find_coordinated_programs` (`R/coordination.R:185+`).
7. **Classes** — `as_cell_program_data` constructor/validator (`R/data-model.R:21-110`: named list of complete numeric obs×genes matrices with row/col names, `sample_metadata` with `observation_id`, no dup ids); `cell_program_fit` is the same list with class `c("cell_program_fit","cell_program_data","list")` appended at `ebmf.R:366` — **there is no separate validator function for the fit class**; downstream code uses `stopifnot(inherits(fit,"cell_program_fit"))` plus ad-hoc checks (`cell-projection.R:87-90`, `ebmf.R:546`, `R/artifact-genes.R:140`).

**Exact seam:** the `run_flash` closure at `R/ebmf.R:290-293` (plus `.fit_via_flashier_utils` at `ebmf.R:467-495` as the template for the argument mapping). A Python-backend adapter is a drop-in `run_flash` that consumes `args` (`data`, `ebnm_fn`, `var_type`, `greedy_Kmax`, `backfit`, `nullcheck`, `verbose`) and returns a flash-like list; everything downstream (fallback ladder, canonicalization, extraction, projection) keys off the returned object's fields.

---

## (d) Required-fields contract for an R `cell_program_fit` built from a Python backend

Per cell type `ct`, `fit$fits[[ct]]` must contain:
- `flash`: list with **`L_pm`** (n_obs × K matrix, rownames = observation IDs) and **`F_pm`** (n_genes × K, rownames = gene IDs) — both NULL only for a skipped cell type; optional `L_psd`, `F_psd`, `F_lfsr` (same shapes; NULL tolerated — see gaps); **`elbo`** (finite scalar; `.fit_ok` `ebmf.R:374` gates the fallback ladder); `n_factors`, `residuals_sd`; `fixed_L`/`fixed_F` only in `covariate_mode="fixed"` (subtracted at `ebmf.R:344`).
- `Y_used` (processed obs×genes matrix), `genes` (character vector), `elapsed_secs`, `fallback` (label), `preproc_replay` = `list(center, scale, residualize_coefficients, design_columns)` (`ebmf.R:256-263`) — required for projection replay (`cell-projection.R:152-155`), `skipped` (TRUE for K=0 sentinel, `ebmf.R:273-281`).
- Top level: `fits`, `preprocessing` (list incl. `feature_sets`, `var_type`, `covariate_mode`, `feature_scale`, `frozen_spec`; `ebmf.R:351-364`), `provenance` (seed/timestamp/versions; `ebmf.R:365-371`), plus inherited `cell_program_data` fields: **`cell_types`**, **`observation_ids`**, **`sample_metadata`** (data.frame with `observation_id`), `pseudobulk`, `feature_metadata`, `shared_genes`, `dims` (`data-model.R:100-109`).
- After `canonicalize_programs()`: `canonicalization` (per-ct `Z`,`W`,`scale_multiplier`,`sign_multiplier`) and **`scores`** (named list, one obs×K matrix per cell type, colnames `<ct>_<k>`) and **`loadings`** (named list, genes×K) — `ebmf.R:602-625`. `project_programs_to_cells` aborts without `scores`/`canonicalization` (`cell-projection.R:88-90`).
- Optional but load-bearing: `attr(flash, "flashbridge_metadata")` — its presence switches canonicalization to the flashier.utils branch (`ebmf.R:559-560`) and is **mandatory for projection** (`cell-projection.R:141-142`). A Python backend that omits it gets the legacy inline canonicalization and loses `project_programs_to_cells` unless a compatible manifest + a projection replacement is supplied.

## (e) Reticulate adapter: what must be translated / what has no equivalent

Must translate:
- **L (loadings/activity, obs×K) and F (factors/effects, genes×K)** posterior means → `L_pm`/`F_pm` with dimnames copied from `Y_used` (rownames/colnames are checked throughout, e.g. `extract.R`, `cell-projection.R`).
- **tau / residual variance**: flashier's `residuals_sd = 1/sqrt(tau)` under `var_type` (scalar, per-row vector, per-col vector, or Kronecker a·b). Python must emit `residuals_sd` in the same shape implied by `var_type`; `given.tau`/`S` (per-gene floor) must be folded in or reported separately.
- **Prior params** (`L_ghat`/`F_ghat`): ebnm `g` objects (e.g. point-laplace a/b mixtures). Python priors (spike-slab, normal, etc.) must be mapped to an ebnm-compatible list or stored as opaque metadata; flashier-utils `prior_identifiers` expects an `ebnm` export name or `"custom_callback_required"` (`fit.R:139-144`).
- **Objective trace**: flashier exposes only final `elbo` on the flash object (per-iteration values are verbose-only); the adapter needs at minimum a finite final ELBO for `.fit_ok` and the fallback ladder; a full trace can be attached as extra metadata but nothing downstream consumes it.
- **K=0 / degenerate sentinel**: return `list(L_pm=NULL, F_pm=NULL, ...)` with `skipped=TRUE` semantics preserved (`ebmf.R:273-281`).

No equivalent / lossy:
1. **lfsr & posterior sd** — `F_lfsr`/`L_psd` come from ebnm posterior summaries. If the Python backend lacks per-element posterior variances/local false sign rates, `top_genes(by="lfsr")` errors (`extract.R:125-128`) and `.cmp_oriented_loading(use_lfsr=TRUE)` silently degrades to unshrunk loadings (`cmp-interpret.R:53-54`); both paths tolerate NULL, so ship NULL and document the degradation.
2. **Native projection oracle** — `project_programs_to_cells` requires `flashbridge_metadata` and calls `flashier.utils::project_factor_activity` on a *native* flash object (`cell-projection.R:141-142,177-186`); a Python fit cannot feed it. Either (i) keep a tiny R-side flashier refit for projection, (ii) implement the projection (regress query on F, `cell-projection.R` replay) in Python and bypass the oracle, or (iii) declare fits unprojectable.
3. **ELBO comparability / fallback semantics** — the fallback ladder and `fit_multistart_ebmf` selection rank candidates by flashier's ELBO (KL terms under ebnm priors); a Python backend's objective is on a different scale, so `.fit_ok` finiteness is portable but *magnitudes are not comparable across backends* — regression tests must not assert cross-backend ELBO equality, and fallback triggers (NaN ELBO) must be re-derived for the Python solver's failure modes.

## (f) Tests: what the ebmf path covers and what a backend-swap test must lock in

`tests/testthat/test-ebmf.R` (123 lines, 6 tests):
1. **End-to-end + canonicalization** (L1-22): planted program in 30×60 data; asserts `scores`/`loadings` shapes & dimnames, and **exact reconstruction preservation**: `max|S W' − L_pm F_pm'| < 1e-8` (L19-21).
2. **Degenerate cell type → K=0 sentinel** (L25-38): `fits$B$skipped`, `ncol(scores$B)==0`, no error.
3. **var_type validation** (L41-47): rejects 3, c(0,1), "1".
4. **Kronecker c(1,2)** (L50-64): fits end-to-end, `preprocessing$var_type == c(1,2)`.
5. **C4 guard** (L67-88): `max_factors > n−q−1` aborts with "residual degrees of freedom", incl. covariate and fixed-mode q.
6. **C4 non-binding + degenerate precedence** (L91-123).

Related coverage: `test-sync-fitting.R` / `test-sync-baseline.R` (flashier.utils backend equivalence), `test-mr1-projection.R`, `test-cmp-interpret.R`, `test-extract-stability.R`, `test-covariates.R`; fixture `helper-toy.R` (`.toy()`: 40×60, two cell types, planted shared program + batch/age effects).

A backend-swap regression test must lock in: (i) reconstruction identity `scores %*% t(loadings) == L_pm %*% t(F_pm)` to ~1e-8 after canonicalization; (ii) dimnames contract (observation IDs on scores rows, gene IDs on loadings rows, `<ct>_<k>` program columns); (iii) K=0 degenerate sentinel behavior; (iv) var_type acceptance/rejection matrix incl. Kronecker; (v) C4 guard firing at `n−q−1`; (vi) finite-ELBO gate + fallback recording (a forced-failure path must land on the recorded `fallback` label, not error); (vii) sign/scale canonicalization determinism (anchor gene, RMS=1) — reproducible under fixed seed; and (viii) if projection is kept: `project_programs_to_cells` replay (residualize→center→scale) and coverage gates. It must NOT assert ELBO equality with the R backend (different objective scale).
