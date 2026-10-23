# cellprograms implementation review

**Repo**: `cellprograms` @ `7c0a212` (branch `main`) — 22 R files, 6,767 lines, 50 exported functions, plus `inst/python/seacells_helpers.py`.
**Date**: 2026-10-19.
**Method**: four parallel module reviews (core fit + data model; projection + sharing + stability; CMP machinery; metacells + plots + extract), every critical/major finding independently re-verified against the code by direct reads before inclusion; minor findings spot-checked for calibration (3/3 spot-checks confirmed). Severity rubric: **critical** = wrong results or crash on valid input; **major** = wrong behavior on edge input, or misleading metadata/docs; **minor** = robustness gap.
**Excluded by construction**: the documented limitation register (L1 projection query-set dependence, L2 CMP C-family FWER 0.040–0.058, L3 sharing n_top/support dilution, L4 frozen spec lacks feature_scale, L5 cross-cohort sign/scale scope, L6 metacells first-use basilisk rebuild, var_floor legacy routing, D-2a forced matches). Those are known and documented; nothing below duplicates them.

## Verdict

The statistical core is sound. Canonicalization math, permutation p-values (`(1+count)/(n_perm+1)` everywhere), blocked permutation, greedy program matching, principal-angle math, the CMP max-T null, seed propagation into Python, ULM/Stouffer/hypergeometric enrichment math, and the absence of `options()` leaks were all verified clean. The findings concentrate in three patterns:

1. **Two genuine correctness bugs** in newer code paths (multi-group soft assignments; covariate-name collision in context association).
2. **Fail-open edge paths** — NA/NaN propagating into `if (NA)` crashes or opaque solver errors — inconsistent with the package's otherwise well-applied `cli_abort` fail-closed convention.
3. **Doc/behavior mismatches** — documented return components, parameter defaults, and claimed refusals that the code does not implement.

Totals: **2 critical, 12 major, 50 minor**.

---

## Critical findings (verified)

### C-1. Multi-group `soft_assignments` is broken — crash or silently wrong matrix
`R/metacells-seacells.R:332-335`

`do.call(rbind, soft_list)` stacks per-group cells×metacells matrices. Metacell columns are group-local: group A's "SEAC_1" and group B's "SEAC_1" are different metacells. When groups produce different metacell counts (the common case — count derives from per-group cell counts), `rbind` errors on mismatched columns. When counts happen to be equal, it silently returns a matrix whose columns mix unrelated metacells across groups (not block-diagonal). Any downstream consumer of `soft_assignments` from a multi-group call gets garbage or a crash. Tests cover only the single-group case, which is why this shipped.

*Fix direction*: build the block-diagonal joint matrix (per-group column blocks, zeros elsewhere) with globally unique metacell names, and add a two-group test with unequal metacell counts.

### C-2. `cmp_context_association` matches covariates named `x*` as context terms
`R/cmp-annotate.R:343`

`rr <- grep("^x", rownames(co))` selects "context" rows of the coefficient matrix by name prefix. A covariate whose name starts with "x" (e.g. `xsite`, `x_batch`) is swept into the context set, so its p-value can be reported as the CMP's context association; because `cmp_triage` takes the best p per CMP, a strong covariate effect can flip a CMP's label. Silent wrong results on valid input.

*Fix direction*: match the exact context term names constructed by the design formula (or an explicit stored term list), never a prefix grep; add a regression test with a covariate named `x*`.

---

## Major findings (verified)

| # | Location | Finding |
|---|----------|---------|
| M-1 | `R/ebmf.R:266-279` | `scale=TRUE` with a constant gene: `scale()` yields a NaN column → `colSums(abs(Y))` is NA → `keep_cols` is NA → `if (!any(keep_cols))` throws "missing value where TRUE/FALSE needed" (or NaN columns propagate to the backend). Valid non-default input crashes opaquely instead of dropping the degenerate gene. Near-critical. |
| M-2 | `R/ebmf_backend.R:117-118,143-147` | cebmf_torch `converged` and `effective_k` are extracted from Python and stored but **never read anywhere in R/**. The only gate is `.cebmf_fit_ok` (finite objective + finite factors), so a non-converged torch fit is silently accepted into results. |
| M-3 | `R/ebmf.R:100` | `@return` documents a `programs` component that is never assigned; `scores`/`loadings` exist only after `canonicalize_programs()`. Users following the docs get `$programs == NULL`. |
| M-4 | `R/cell-projection.R:117` | `fs <- if (!is.null(feature_scale)) feature_scale else fit$preprocessing$feature_scale` — a user-supplied `feature_scale` is adopted without comparing to the fit's recorded units, although the roxygen claims a mismatched query "is refused". Silent wrong-space projection when the user passes a mismatched value. Near-critical. |
| M-5 | `R/stability.R:74-75` | `sample.int(nrow(Y), max(5L, floor(resample_frac*nrow(Y))))` requests ≥5 rows; a cell type with <5 donors crashes `sample.int`. This runs outside the CP-1 tryCatch failure ledger, so one small cell type kills the whole stability run unrecorded. |
| M-6 | `R/plots.R:82` | `plot_program_scores` does `do.call(cbind, ...)` over per-cell-type score matrices, requiring equal row counts; the data model explicitly permits per-type observation sets, so this crashes on valid fits. |
| M-7 | `R/plots.R:347` | `plot_sharing_spectrum` uses `scale_color_brewer(palette = "Set1")` (max 9 colors) with one label per pair; >9 distinct pairs crashes. |
| M-8 | `R/plots.R:368-378` | `plot_artifact_flags` hardcodes `vdj_frac`/`genotype_marker_frac` columns; output of `flag_artifact_genes(fit, sets = "sex")` (a documented option) crashes it. |
| M-9 | `R/extract.R:195-201` | `union_panel_loadings`: `scale()` on a constant gene yields NaN → `lm.fit` aborts "NA/NaN/Inf in 'y'"; the `G[is.nan(G)] <- NA_real_` cleanup sits after the abort and is unreachable. |
| M-10 | `R/artifact-genes.R:53` | VDJ regex `"^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD][VC])"` matches **no J-segment genes** (IGHJ/IGKJ/IGLJ/TRAJ/TRBJ/TRGJ/TRDJ all unmatched) and misses IGLL1/IGLL5, despite the doc claiming V(D)J coverage. Under-flagging of the most common B/T-cell artifact class. |
| M-11 | `R/metacells-seacells.R:126,187`; `inst/python/seacells_helpers.py:101-103` | Doc/behavior mismatches: roxygen says `max_iter` default 100 but the signature and Python use 20; `use_sparse` is threaded through R but silently dropped by the Python helper (SEACells 0.3.3 has no such kwarg). |
| M-12 | `R/metacells-seacells.R:145-147` vs `354-360` | `@return` promises mcRigor QC columns in the metacell object's metadata; the constructed object carries only `Group`, `n_cells`, and the `group.by` columns. QC lives only in `out$qc`. |

### Downgraded during verification

- **Fixed-mode `preproc_replay` incompleteness** (reported as major): **not a bug**. `project_programs_to_cells()` refuses fixed-covariate and legacy fits loudly with a clear `cli_abort`, and the refusal **is** documented in that function's roxygen (`cell-projection.R:23-25,154-157`). Residual minor gaps (counted below): `fit_programs()`'s `covariate_mode`/`backend` docs don't warn that fixed-mode and cebmf_torch fits are unprojectable, and the abort message ("legacy or fixed-covariate backend") mislabels cebmf_torch fits, which also lack the flashbridge manifest.

---

## Minor findings (50)

Spot-checked for calibration (validate 1-observation crash, `min_donors` row-count, flip=0 sign concordance — all confirmed); the remainder are subagent-reported with file:line evidence, not individually re-verified.

**Core fit + data model**
1. `identical()` type-strict frozen-spec conflict check can miss integer/numeric mismatches (`ebmf.R`).
2. `n_variable_genes` not validated (negative/NA accepted until downstream failure) (`ebmf.R`).
3. `max_factors` not validated (`ebmf.R`).
4. Undocumented ≥10 residual-df abort in the covariate path (`ebmf.R`).
5. `validate_cell_program_data` crashes on a 1-observation cell type: `var` of length-1 columns → NA → `if (NA)` error, pre-empting its own small-N warning (`data-model.R:129-135`). **[verified]**
6. `feature_metadata` accepted unvalidated.
7. Sign-anchor tie-break is locale-dependent (gene-name sort).

**Projection**
8. No zero-scale guard in the replay `sweep(..., "/")` (`cell-projection.R`).
9. `min_coverage`/`maxiter` unchecked.
10. Fit-side docs silent on unprojectability of fixed-mode/cebmf fits; abort message mislabels cebmf fits (see "Downgraded" above). **[verified]**

**Sharing / principal angles / stability / coordination**
11. `set.seed` leaks (global RNG mutation) at `program-sharing.R:111`, `principal-angles.R:72`, `stability.R:55`, `coordination.R:319`, `cmp-interpret.R:149`.
12. Degenerate-program NA propagation in sharing scores.
13. `.sharing_matrices` called twice per `program_sharing()` run (wasted compute).
14. 0-shared-observation pair → opaque SVD failure (`principal-angles.R:169-175`).
15. `sharing_spectrum` default unseeded.
16. `.pa_engine` `na.rm` handling biases the null (`:85`).
17. `match_programs` space documentation gap.
18. `conditional_recovery_freq` dilution (`stability.R:199`).
19. Legacy fits: `provenance$seed` is `numeric(0)` (`:101`).
20. `min_donors` filter counts matrix rows, not observed donors — NA rows inflate the count (`coordination.R:79`). **[verified]**

**CMP machinery**
21. Sign-concordance edge: tie sum (`flip = 0`) deterministically scored discordant; an all-NA sign vector would score concordant = TRUE (`cmp-stability.R:236-237`). **[verified]**
22. Unnamed `cmp_args` bind positionally (`cmp-stability.R:120`).
23. Dead `set.seed` (`cmp-stability.R:143`).
24. `ambient_frac` NA → 0 (`cmp-interpret.R:341`).
25. Hardcoded IEG threshold (`cmp-interpret.R:404`).
26. NA `mor` propagation (`cmp-interpret.R:163`).
27. `prepare_collectri` integer 0/1 → all `+1` mor (`cmp-annotate.R:44-46`).
28. `cmp_gene_overlap` returns NULL on empty input (`:347`).
29. `donor_check` output columns conditional on input (`:151`).

**Metacells**
30. No warning listing skipped groups (`metacells-seacells.R:113-114`).
31. All-zero soft-assignment rows map to `mc_names[1]` (`:408`).
32. `tmp_root` collision risk across concurrent runs (`:236`).
33. Metadata column collisions when a `group.by` name shadows `Group`/`n_cells` (`:356-360`).
34. Fallback QC column mismatch (`:296-299`).
35. mcRigor availability guard swallowed by tryCatch (`metacells-mcrigor.R:229-231`).
36. NA verdict coerced to `dubious = TRUE` (`:296-299`).
37. `.mcrigor_optimize_pick` index misalignment (`:318-326`).
38. Membership doc type mismatch (`:236-239`).
39. Silent group drops (`:213-216`).

**Metacells IO**
40. Layer fallback order undocumented (`metacells-io.R:9-19`).
41. `dims` unvalidated (`:22-27`).
42. Numeric cell names → Python KeyError (`:44-47`).

**Plots**
43. `se = Inf` at t = 0 (`plots.R:740`).
44. `plot_variance_explained` positional PVE mapping (`:209-214`).
45. Annotation df without `observation_id` → all-NA join, silent (`:60-63`).

**Extract**
46. Rank-deficient fixed-mode replay → `solve()` crash (`extract.R:186-188`).
47. Zero-overlap gene panel unhandled (`:169-172`).

**Mofacellular**
48. Name ambiguity via `paste0(ct, "_", donor)` when names contain `_` (`mofacellular.R:41-44`).
49. Duplicate `cell_counts` keys (`:70-73`).
50. 0-row pseudobulk unhandled (`:29-31`).

---

## Verified clean (no action)

- Canonicalization math (RMS scaling + anchor-feature sign) exact; matches `flashier.utils::standardize_factors` to 0 difference.
- Permutation p-values use `(1+count)/(n_perm+1)` at every site; blocked permutation permutes within strata on the ct_b side only; CMP max-T null permutes coherently within cell type.
- Greedy program matching index-correct; principal-angle math correct.
- Seed propagation into Python SEACells correct; no `options()` leaks.
- ULM / Stouffer / hypergeometric enrichment math correct.
- CMP statistics core overall reviewed clean (consistent with the 5500-replicate calibration evidence).

## Test-coverage gaps exposed by the findings

- Multi-group `use_soft_assignments` untested (C-1 shipped because of this).
- No test with a covariate named `x*` (C-2).
- No test with `scale=TRUE` + constant gene (M-1), `<5`-donor stability (M-5), `sets="sex"` plotting (M-8), >9-pair sharing spectrum (M-7).
- The opt-in `CPMC_INTEGRATION=1` metacells end-to-end test has never been run in any suite.

## Recommended fix priority

- **P0 (correctness)**: C-1 block-diagonal soft assignments; C-2 exact-term matching in `cmp_context_association`.
- **P1 (crash on valid/edge input)**: M-1 NaN-tolerant `keep_cols`; M-5 wrap resampling in the failure ledger + guard n<5; M-6/M-7/M-8 plotting guards; M-9 drop constant genes before `lm.fit`.
- **P2 (silent quality/metadata)**: M-2 gate or warn on `converged == FALSE`; M-4 enforce or drop the documented refusal; M-10 extend the VDJ regex to J segments + IGLL1/5; M-3/M-11/M-12 doc fixes.
- **P3**: the 50 minors, batched by module.

Suite state at review time: FAIL 0 | WARN 3 | SKIP 1 (cellprograms), FAIL 0 | WARN 11 | SKIP 6 (flashier-utils) — all findings above are invisible to the current suites, which is itself the main process takeaway: the fail-closed convention and the test matrix need to extend to the newer modules.
