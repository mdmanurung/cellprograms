# Regression guard for .match_programs_loadings(): p-values sit on a 1/(n_perm + 1) grid and
# depend on every permuted correlation, so identical p pins the RNG stream and the null
# distribution. Golden captured from the cor()-on-permuted-copies loop (release/0.3.0 @ 0914144).
# stats::cor() forces `y` before `x`, so that loop draws Wb's permutation first. Loadings with
# exact ties are deliberately not pinned: there p depends on the platform's floating-point
# arithmetic (macOS CI broke ties differently), and a crossprod of pre-scaled columns changed p
# in 36 of 180 fuzzed configurations; see benchmarks/performance/RESULTS.md.

golden <- readRDS(test_path("fixtures", "sharing_perf_golden.rds"))

test_that(".match_programs_loadings output is unchanged", {
  got <- .match_programs_loadings(.sharing_perf_fit(), n_perm = 199L, seed = 1L)
  want <- golden$continuous
  expect_identical(got$p, want$p)
  expect_identical(got[c("ct_a", "prog_a", "ct_b", "prog_b")], want[c("ct_a", "prog_a", "ct_b", "prog_b")])
  expect_equal(got$cor, want$cor, tolerance = 1e-10)
  expect_equal(got$q, want$q, tolerance = 1e-12)
  # the fixture is not degenerate: shared program found, p-values spread
  expect_equal(want$p[1], 1 / 200)
  expect_gt(length(unique(want$p)), 20L)
})
