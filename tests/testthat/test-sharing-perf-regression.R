# Regression guard for .match_programs_loadings(): p-values sit on a 1/(n_perm + 1) grid and
# depend on every permuted correlation, so identical p pins the RNG stream and the null
# distribution. Golden captured from the cor()-on-permuted-copies loop (release/0.3.0 @ 0914144).
# stats::cor() forces `y` before `x`, so that loop draws Wb's permutation first; the second
# fixture has discrete loadings, where exact ties make p depend on the arithmetic used for
# the null (a crossprod of pre-scaled columns changed p in 36 of 180 fuzzed configurations).

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

test_that("p-values with exactly tied correlations are unchanged", {
  got <- .match_programs_loadings(.sharing_tie_fit(), n_perm = 99L, seed = 1L)
  expect_identical(got$p, golden$ties$p)
  expect_identical(got$q, golden$ties$q)
  expect_gt(length(unique(golden$ties$p)), 2L)
})
