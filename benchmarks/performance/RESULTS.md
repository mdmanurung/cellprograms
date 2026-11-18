# Performance of the permutation loops

Measured with `bench_cmp.R` (best of 3 wall-clock seconds; the script also times `.match_programs_loadings()`, 3.0 s for 10 cell types, 2000 genes, `n_perm = 199`, unchanged), R 4.5.1, one core, synthetic
planted scores. Before = `0914144` (release/0.3.0), after = the `perf/cmp-fast-path` branch.
These are timings on generated data, not a result about biology.

| Call | Size | Before | After | Ratio |
|---|---|---|---|---|
| `find_coordinated_programs()` | 5 cell types, 200 donors, `n_perm = 99` | 1.06 s | 0.24 s | 4.4x |
| `find_coordinated_programs()` | 10 cell types | 2.31 s | 0.84 s | 2.8x |
| `find_coordinated_programs()` | 20 cell types | 10.49 s | 3.45 s | 3.0x |

A single-run profile at 10 cell types, 100 donors, `n_perm = 999` gave 25.7 s before and
5.4 s after.

## What changed

- `find_coordinated_programs()`: the loop rebuilt a labelled data frame per cell-type
  pair and permutation and used only its `stat` column (`.cmp_edges()` was 97% of the run
  time). It now uses a precomputed plan, a numeric-only edge statistic with the same
  arithmetic, and an integer-indexed node maximum. Results are identical, not close:
  `node_p`, edge and node statistics match the old code on every pinned case.

## Tried and not adopted: loading match

Scaling the loading columns once and taking one `crossprod()` per permutation made
`.match_programs_loadings()` 1.5x faster (3.0 s to 2.0 s) and gave identical p-values on
continuous loadings. A differential fuzz run found that it changes p-values when loadings
contain exact ties (discrete values; 36 of 180 configurations, p off by up to 0.38),
because `cor()` and a cross-product of scaled columns differ by about 1e-17, which breaks
ties in `|null| >= |obs|` differently. Real EBMF loadings are continuous, but 1.5x is not
worth a change in seeded results, so the code is unchanged. The golden test on
continuous loadings stays. A tie fixture was tried and dropped: on exact ties p depends on the
platform's floating-point arithmetic (macOS CI broke ties differently from Linux), so it cannot
be pinned across platforms. The fuzz run also showed that `stats::cor()` forces `y` before `x`, so
the old loop draws the second matrix's permutation first; any rewrite must keep that order.

## What is left, and why there is no Rcpp

- `find_coordinated_programs()`: `stats::cor(use = "pairwise.complete.obs")` is 36% of
  the remaining time, and R already runs it in C. The rest is R-level glue per pair.
  Rcpp could cut the 36% by at most about half, so a compiled dependency would buy
  roughly 20% for a harder install. The plan's gate (go to Rcpp only if `cor` plus the
  group maximum is half the time or more) says stop.
- `.match_programs_loadings()`: with the crossprod change, 54% of the time was `sample.int()`; the seeded stream pins those draws.
- Not worth work (measured): `ambient_program_scores()` is 1.4% of a
  `find_coordinated_programs()` call and 1.0-1.7% of `assess_cmp_stability()`;
  `annotate_cmp_lr()` takes about 1 s.
- `assess_program_stability()` and the other fit-based functions spend their time in
  `flashier`/`ebnm` (212-870 s per fit in `benchmarks/coordination/PREREG.md`); R-level
  changes cannot help there.

`assess_cmp_stability()` calls `find_coordinated_programs()` once per replicate, so it gets
the same factor. Not measured separately.

## Regression guards

`tests/testthat/test-cmp-perf-regression.R` and `test-sharing-perf-regression.R` pin the
outputs captured at `0914144` (`tests/testthat/fixtures/`, with the scripts that made them).
Do not regenerate them to make a failing test pass.
