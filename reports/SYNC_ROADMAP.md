# Cross-Repository Synchronization — Stage A Receipt

**Status:** Stage A (SYNC-0 → SYNC-3) implemented on the `sync/stage-a` branches of
both repositories. This document records the supported version pairs, the parity
evidence, known deviations, and the rollback plan. It links from, and does not
replace, the scientific `reports/ROADMAP.md`.

## Dependency direction

```
flashier + ebnm  →  flashier.utils  →  cellprograms
```

One-way. Scientific EBMF operations (fit, representation, projection, uncertainty)
live upstream in `flashier.utils`; cell-type orchestration, donor structure, CMPs
and biology live downstream in `cellprograms`. Never the reverse.

## Supported version pair (Stage A)

| Repo | Branch | Commit | Key version |
|---|---|---|---|
| `flashier` (dep) | pinned | `fd44811ab22d61ee905e1964175571c452fe7683` | 1.0.60 (dev) |
| `flashier-utils` | `sync/stage-a` | SYNC-1A commit | 1.1.0.9000 |
| `cellprograms` | `sync/stage-a` | SYNC-2 commit | 0.2.0 |

**Critical environment note (SYNC-0 finding):** `flashier.utils` is sensitive to the
exact `flashier` commit. Against CRAN flashier 1.0.7 the upstream suite showed 16
failures + 5 erroring files (backfit `get.obj` errors, missing `psych` engine);
against the pinned `flashier@fd44811` (1.0.60) + `psych`, the suite is clean
(1132 passed / 0 failed). The `Remotes` pin is therefore load-bearing, not cosmetic.

## What changed

### SYNC-0 — freeze contracts
- `flashier-utils`: added signed, tied-anchor, zero-activity, zero-factor fixtures
  (`tests/testthat/helper-fixtures.R`). Baseline 1114 → 1132 passed / 0 failed.
- `cellprograms`: added `helper-sync-fixtures.R` + `test-sync-baseline.R` freezing the
  ordinary-fit output contract (fields, exact reconstruction, unit-RMS convention,
  max-|loading| sign anchor, `<cell_type>_<k>` IDs, seeded fingerprint). 32 assertions.

### SYNC-1A — opt-in RMS upstream (`flashier-utils`)
- `standardize_factors(scaling = "rms")`: unit-RMS activity, compensating effects,
  reconstruction preserved exactly; `s_k = sqrt(mean(A[,k]^2))`, zero-activity guard
  `s_k = 1`; per-factor scale recorded as `manifest$rms_scale`. Defaults (`raw`/`ldf`)
  unchanged. 18 parity assertions in `test-representation-rms.R`.

### SYNC-1B — ordinary-fit migration (`cellprograms`)
- Ordinary fits (no fixed covariates) route through `flashier.utils::fit_ebmf()` with an
  exact argument mapping. **Verified to machine precision** against the legacy
  `do.call(flashier::flash, ...)` path: identical `n_factors`, ELBO, `L_pm`, `F_pm`, and
  canonicalized scores/loadings (max abs diff 0).
- Fallback ladder (point_laplace → point_normal; nullcheck-off), per-attempt seed reset,
  non-finite-ELBO rejection, and `$fallback` recording unchanged.
- Unmappable `flash_control` overrides are **rejected explicitly** (error naming the
  setting and the opt-out), never silently dropped.
- `.flash_fixed()` (fixed-covariate path) and degenerate K=0 sentinels stay on the legacy
  backend, unchanged. `options(cellprograms.use_flashier_utils = FALSE)` forces legacy.
- 27 assertions in `test-sync-fitting.R`.

### SYNC-2 — canonicalization parity (`cellprograms`)
- `canonicalize_programs()` routes ordinary fits through
  `standardize_factors(scaling = "rms", orientation = "anchor_feature")`, reproducing the
  inline RMS + sign convention to 0 difference; native dimnames preserved on Z/W.
- Native `F_psd` / `F_lfsr` / `scale_multiplier` / `sign_multiplier` semantics preserved;
  canonicalization idempotent; program-sharing `F_psd * scale_multiplier` invariant holds.
- 22 assertions in `test-sync-representation.R`.

### SYNC-3 — CI, version compatibility, release
- `cellprograms/DESCRIPTION`: R ≥ 4.3; `flashier`, `ebnm`, `flashier.utils` moved to
  `Imports`; `Remotes` pins flashier + flashier-utils. R CMD check: **Status OK**.
- `cellprograms` CI: pinned-compatible lane (ubuntu + macOS) plus an `upstream-main` lane
  tracking `flashier-utils@main` to catch upstream drift.
- `flashier-utils` CI: new `revdep-cellprograms.yaml` runs the downstream sync/core tests
  on upstream PRs and weekly.

## Parity evidence (the release gate)

| Check | Result |
|---|---|
| `fit_ebmf()` vs direct `flash()` (toy, K=3) | identical ELBO, L_pm, F_pm (diff 0) |
| Upstream RMS vs legacy canonicalization | scores/loadings diff 0 |
| cellprograms full suite | **438 passed / 0 failed** |
| flashier-utils full suite | **1132 passed / 0 failed** |
| cellprograms R CMD check (`--no-manual`) | **Status OK** |

## Known deviations / unsupported settings

- `flash_control` members beyond the core set (`data, ebnm_fn, var_type, greedy_Kmax,
  backfit, nullcheck, verbose`) are not mapped through `fit_ebmf()`; they error with
  instructions to use the legacy backend. This is deliberate (no silent control loss).
- The fixed-covariate path (`covariate_mode = "fixed"`, `.flash_fixed()`) is **not**
  migrated and remains on the legacy backend by design.
- R ≥ 4.3 is now required (was ≥ 4.1), following `flashier.utils`. R 4.1 is no longer
  supported once the hard dependency is in `Imports`.

## Rollback plan

- Set `options(cellprograms.use_flashier_utils = FALSE)` to force the legacy ordinary-fit
  backend at runtime (no reinstall needed).
- The legacy code paths (direct `flashier::flash`, inline RMS canonicalization) are
  retained in `R/ebmf.R` and exercised by the test suite, so reverting the migration is a
  matter of changing the default of that option or pinning `cellprograms` to the commit
  before SYNC-1B.
- `flashier.utils` `raw`/`ldf` scalings are untouched; the `rms` addition is purely
  additive and can be ignored by existing callers.

## Not in scope (Stage B)

Multiresolution discovery (MR-0 → MR-4: frozen preprocessing, PB→cell projection, cNMF
baseline, metacells, mcRigor, anchored refinement, ablation benchmarking) and the OPT
utility consolidation are **not** part of Stage A and are gated behind this release.
