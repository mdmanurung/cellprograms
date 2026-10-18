# Covariates Audit — `cellprograms` (final report)

**Scope**: covariate functionality of `cellprograms` 0.2.0 — `.cov_design`, `.flash_fixed`,
`fit_celltype_programs` covariate paths (residualize / fixed), `preprocessing_spec` interplay,
`project_programs_to_cells` replay, stability refits, `union_panel_loadings`.
**Method**: adversarial test suite `/workspace/cov_audit_tests/test_covariates.R` (46 checks),
`trace()` instrumentation of `.flash_fixed`, targeted algebraic reconstruction of the fit's
internal matrices. Repo `/workspace/repo2_full` was **not modified** (verified: git clean,
line numbers unchanged from first pass).

**Final suite result: 46 PASS, 0 FAIL** (exit 0). Failures of defect-marker tests are
intentional (they pass *while* the defect is present; labels say "invert when fixed").

## 1. Findings table

| ID | Class | Severity | Code | Evidence (test) |
|---|---|---|---|---|
| F1 Replay center/scale misaligned when zero-variance genes are dropped | CONFIRMED DEFECT | **P2** | `R/ebmf.R:211-219` (replay set from full panel), `R/ebmf.R:246-249` (drop block does not subset replay), `R/cell-projection.R:174-176` (projection indexes by post-drop positions) | T17a: `length(center)=150` vs `ngenes=149`; T17b: code-replayed centering off by **max 10.43** vs `Y_used`, correctly-aligned replay exact (0.00e+00). Affects the default (no-covariate) path too. |
| F2 Rank-deficient (collinear) design poisons projection replay | CONFIRMED DEFECT | **P2** | `R/ebmf.R` residualize block keeps aliased columns (no pivot); `R/cell-projection.R:168-170` (`X %*% beta` with NA coefficients) | T16a: `age_copy` row of `residualize_coefficients` is NA; T16b: projection errors "Missing values are not supported for this operation" (fail-closed, but a legitimate collinear design is unusable at projection time) |
| F3 Frozen spec bypasses fixed-mode centering invariant | CONFIRMED DEFECT | **P2** | guard at `R/ebmf.R:116-117` runs *before* spec substitution (`R/ebmf.R:157-181`) | T14f: a spec frozen with `covariate_mode="fixed", center=FALSE` refits without error — fixed mode ran with `center=FALSE` |
| F4 Interaction term referencing a within-ct-constant covariate crashes `model.matrix` obscurely | CONFIRMED DEFECT | P3 | `R/ebmf.R:398-402` (drops the variable from `m` but leaves `age:sex` in the formula) | T21: error "variable lengths differ (found for 'sex')" — obscure, not actionable |
| F5 Obscure (but fail-closed) projection errors on degenerate query metadata | DOCUMENTED LIMITATION | P3 | `R/cell-projection.R` `.replay_design` | T13b: lost covariate level aborts with model.matrix's "contrasts can be applied only to factors with 2 or more levels" (not the design-mismatch message); T13e: missing column → "undefined columns selected". Safety property (abort) holds in both. |
| F6 Adjustment ordering: variable-gene selection on RAW Y; residualize → center → scale; `scale=TRUE` records residual SDs | VERIFIED-CORRECT | — | `R/ebmf.R:203-219` | T1a/T1b/T1c (diffs 0.00e+00; scale−residSD 4.4e-16 vs rawSD 2.10) |
| F7 df guard uses design **rank** − 1 | VERIFIED-CORRECT | — | `R/ebmf.R:255-273` | T2a/T2b: collinear `age+age_copy` → q=1; max_factors=10 fits, 11 trips guard |
| F8 Residualization orthogonal to design; lm.fit-based; coefficients stored | VERIFIED-CORRECT | — | `R/ebmf.R:221-224` | T3: max\|t(R)X\| = 7.1e-15 (full-rank), 1.8e-14 (collinear) |
| F9 `.cov_design` semantics: formula ≡ character, droplevels, constant-within-ct dropped/not-adjusted, NA abort, nrow−rank≥10 guard, leverage warning, fixed-mode intercept drop + center/scale + rank pivot | VERIFIED-CORRECT | — | `R/ebmf.R:377-427` | T4, T5, T6a/b/c, T7, T8, T9, T10a/b, T14e, T15, T22 |
| F10 `.flash_fixed` semantics: `fixed_L` = design exactly; `fixed_F` = OLS init (backfit=FALSE); Kmax counts **added** factors so `max_factors` = free programs in both modes; fixed columns not returned as free programs; planted covariate module removed (residualize and fixed); fixed-mode fits refuse projection | VERIFIED-CORRECT | — | `R/ebmf.R:429-464`, `R/cell-projection.R` manifest check | T11a (identical), T11b (0.00e+00, §2), T11c (K_free=2, names B_1,B_2), T11d (\|cor(age,score)\|: 0.960 → 0.001 both modes), T11e, T18; `cor(fixed_F age col, planted loadings)=0.991` |
| F11 Projection replay reproduces in-sample adjustment exactly; order-independent; design-mismatch fail-closed | VERIFIED-CORRECT | — | `R/cell-projection.R:150-180` | T12a/b (3.6e-15), T12d, T13a/c/d/e |
| F12 Frozen-spec conflict detection (explicit args refused; digest recorded; rehydrated formula) | VERIFIED-CORRECT | — | `R/ebmf.R:135-181` | T14a–e |
| F13 Stability refits and `union_panel_loadings` carry covariates | VERIFIED-CORRECT | — | `R/stability.R`, `R/program-sharing.R` | T19 (8 programs), T20 (no errors, finite) |
| F14 T11b "fixed-mode Y_used/fixed_F inconsistency (6.54)" | FALSE-POSITIVE (harness bug) | — | test used `scale(Yraw)` — base `scale()` defaults to `scale=TRUE`; the fit uses `scale=FALSE` | §2 below |

## 2. T11b resolution — FALSE ALARM (test-harness bug, not a package defect)

**Outcome: the fixed-mode pipeline is exactly self-consistent. No P1.**

Resolved by `trace()`ing `cellprograms:::.flash_fixed` (dump list initialized in a
`globalenv()` environment *before* `trace()`, per the delegated method) plus algebraic
reconstruction:

1. `.flash_fixed` receives **exactly** `scale(panel, center=TRUE, scale=FALSE)` as `args$data`
   (max|diff| = 0.000e+00; cell-type order B,CD4,CD8,NK,Mono; first call = B).
2. `Y_used == scale(panel, center=TRUE, scale=FALSE) − fixed_L %*% t(fixed_F)` to
   **6.7e-16** (`R/ebmf.R:334` subtracts the fixed part before storing `Y_used`).
3. `fixed_F` **is** the OLS initializer exactly for `backfit=FALSE`
   (max|fixed_F − OLS| = 0.000 against the *center-only* panel). The first pass's claim is
   confirmed; my interim diagnostic value 1.912 was an artifact of comparing against the
   wrongly-scaled panel (same bug as T11b itself).
4. Root cause of the 6.54: the test computed `Yc <- scale(Yraw)` — base `scale()` defaults to
   `scale=TRUE` — and compared a center+**scaled** panel against a center-only pipeline.
   max|center+scale − center-only| = 6.538 = the anomaly, bit-for-bit. Fits are bit-identical
   across processes (sum(fixed_F) = 5.319621 in every run), ruling out nondeterminism.

The corrected T11b (and a note line documenting the fixed_F/OLS relationship) is in the
final test script; it passes at 0.00e+00.

## 3. Proposed patches (unified diffs — NOT applied to the repo)

### F1 / P2-1 — subset replay center/scale in the gene-drop block (`R/ebmf.R:246-249`)

```diff
     if (!all(keep_cols)) {
       cli::cli_inform("Dropping {sum(!keep_cols)} zero-variance gene(s) in cell type {.val {ct}}.")
       Y <- Y[, keep_cols, drop = FALSE]
       genes <- genes[keep_cols]
+      replay$center <- replay$center[keep_cols]
+      replay$scale  <- replay$scale[keep_cols]
     }
```

### F2 / P2-2 — zero NA coefficients before projection replay (`R/cell-projection.R:168-170`)

```diff
       beta <- replay$residualize_coefficients
+      # Aliased columns of a rank-deficient design carry NA coefficients (lm.fit pivot);
+      # the in-sample fitted values used the pivoted solution, i.e. NA <=> 0.
+      beta[is.na(beta)] <- 0
       Q <- as.matrix(expression[cells_ct, shared, drop = FALSE]) -
         X %*% beta[replay$design_columns, shared, drop = FALSE]
```

### F3 / P2-3 — repeat the fixed-mode centering guard after spec substitution (`R/ebmf.R`, insert after line 181)

```diff
     frozen_digest <- preprocessing_spec$digest
   }
+
+  # The guard above ran before spec substitution; a spec frozen with
+  # covariate_mode = "fixed", center = FALSE must not bypass it.
+  if (covariate_mode == "fixed" && !is.null(covariates) && !center) {
+    cli::cli_abort("{.code covariate_mode = \"fixed\"} requires {.code center = TRUE}.")
+  }
```

### F4 / P3 — drop interaction terms that reference a within-ct-constant covariate (`R/ebmf.R:398-402`)

```diff
   if (any(const)) {
     cli::cli_inform("Covariate(s) {.val {vars[const]}} constant in cell type {.val {ct}}; dropped for this cell type.")
-    f <- stats::drop.terms(stats::terms(f), which(attr(stats::terms(f), "term.labels") %in% vars[const]),
-                           keep.response = FALSE)
+    tl <- attr(stats::terms(f), "term.labels")
+    tl <- tl[vapply(tl, function(t) !any(all.vars(stats::as.formula(paste("~", t))) %in% vars[const]),
+                    logical(1))]
+    f <- stats::reformulate(tl)   # all-constant case already returned NULL above
     m <- m[, vars[!const], drop = FALSE]
   }
```

## 4. Test-harness bugs fixed this pass (all were harness, not package)

- **T5**: compared against a subset of the full design with `identical()` — fails on
  `assign`/`contrasts` attributes; now compared (attribute-insensitively) against an
  independent `droplevels`ed sex-only `model.matrix`.
- **T6c**: cli wraps the message between "not" and "adjusted"; pattern is now `not\s+adjusted`.
- **T11b**: `scale(Yraw)` defaults to `scale=TRUE` — see §2.
- **T11c/d**: `$scores` is NULL before `canonicalize_programs()`; canonicalized fits are now
  created once in T11 and reused by T12/T18 (re-canonicalizing the same object was also
  removed).
- **T12c/d**: the projection's `activity` matrix intentionally spans **all** cell types'
  programs (8 cols here); comparisons now use B's named program columns.
- **T13b**: assertion relaxed to the safety property (fail-closed abort); the obscure
  contrasts message is recorded as F5.
- **T20**: `fin()` used scalar `||` on matrices; now vectorized `|`, list-safe.
- **T14f/T21**: converted to explicit defect markers (pass while F3/F4 present).

## 5. Overall verdict on scientific soundness

**Sound.** The covariate machinery does what it claims, verifiably and to machine precision:
adjustment is applied in the documented order on the documented panel; residuals are exactly
orthogonal to the design; the degrees-of-freedom guard is rank-aware; the fixed-covariate mode
pins the design as fixed loadings with the OLS solution and returns only free programs;
projection replay reproduces the in-sample adjustment exactly and fails closed on every
metadata mismatch tested. No P0/P1 defect exists; nothing silently corrupts a well-specified
fit. The three P2 defects are edge-case robustness gaps — (F1) silently misaligned replay when
the panel contains zero-variance genes (the most serious: wrong projected activities, default
path included), (F2) collinear designs unusable at projection time (fail-closed), (F3) a
guard-ordering hole reachable only via a deliberately frozen bad spec — plus one P3 crash and
two P3 message-quality notes. All four proposed patches are small, local, and behavior-
preserving elsewhere; applying them before release is recommended.

## Files

- `/mnt/results/covariate_audit/report.md` — this report
- `/mnt/results/covariate_audit/test_covariates.R` — final audit suite (46 checks)
- `/mnt/results/covariate_audit/test_output.txt` — complete suite output (46 PASS, 0 FAIL)
