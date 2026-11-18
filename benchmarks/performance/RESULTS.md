# Performance of the permutation loops

Measured with `bench_cmp.R` (best of 3 wall-clock seconds), R 4.5.1, one core, synthetic
planted scores. Before = `0914144` (release/0.3.0), after = the `perf/cmp-fast-path` branch.
These are timings on generated data, not a result about biology.

| Call | Size | Before | After | Ratio |
|---|---|---|---|---|
| `find_coordinated_programs()` | 5 cell types, 200 donors, `n_perm = 99` | 1.06 s | 0.24 s | 4.4x |
| `find_coordinated_programs()` | 10 cell types | 2.31 s | 0.84 s | 2.8x |
| `find_coordinated_programs()` | 20 cell types | 10.49 s | 3.45 s | 3.0x |
| `.match_programs_loadings()` (inside `program_sharing()`) | 10 cell types, 2000 genes, `n_perm = 199` | 3.02 s | 2.04 s | 1.5x |

A single-run profile at 10 cell types, 100 donors, `n_perm = 999` gave 25.7 s before and
5.4 s after.

## What changed

- `find_coordinated_programs()`: the loop rebuilt a labelled data frame per cell-type
  pair and permutation and used only its `stat` column (`.cmp_edges()` was 97% of the run
  time). It now uses a precomputed plan, a numeric-only edge statistic with the same
  arithmetic, and an integer-indexed node maximum. Results are identical, not close:
  `node_p`, edge and node statistics match the old code on every pinned case.
- `.match_programs_loadings()`: scale the loading columns once, then one `crossprod()`
  per permutation. Seeded p-values are identical.

## What is left, and why there is no Rcpp

- `find_coordinated_programs()`: `stats::cor(use = "pairwise.complete.obs")` is 36% of
  the remaining time, and R already runs it in C. The rest is R-level glue per pair.
  Rcpp could cut the 36% by at most about half, so a compiled dependency would buy
  roughly 20% for a harder install. The plan's gate (go to Rcpp only if `cor` plus the
  group maximum is half the time or more) says stop.
- `.match_programs_loadings()`: 54% of the time is `sample.int()`. The seeded stream must
  not change, so those draws stay.
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
