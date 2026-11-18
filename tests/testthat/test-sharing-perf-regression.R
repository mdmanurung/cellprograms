# Regression guard for .match_programs_loadings(): p-values sit on a 1/(n_perm + 1)
# grid and depend on every permuted correlation, so identical p pins the RNG stream
# and the null distribution. Golden captured before the loop was rewritten.

test_that(".match_programs_loadings output is unchanged", {
  golden <- readRDS(test_path("fixtures", "sharing_perf_golden.rds"))
  got <- .match_programs_loadings(.sharing_perf_fit(), n_perm = 199L, seed = 1L)
  expect_identical(got$p, golden$p)
  expect_identical(got[c("ct_a", "prog_a", "ct_b", "prog_b")], golden[c("ct_a", "prog_a", "ct_b", "prog_b")])
  expect_equal(got$cor, golden$cor, tolerance = 1e-10)
  expect_equal(got$q, golden$q, tolerance = 1e-12)
  # the fixture is not degenerate: shared program found, p-values spread
  expect_equal(golden$p[1], 1 / 200)
  expect_gt(length(unique(golden$p)), 20L)
})

test_that("the permuted correlation equals cor() of row-permuted copies", {
  set.seed(2)
  G <- 120
  Wa <- matrix(rnorm(G * 3), G, 3)
  Wb <- matrix(rnorm(G * 4), G, 4)
  set.seed(9); pa <- sample.int(G); pb <- sample.int(G)
  expect_equal(.perm_cor(.scale_cols(Wa), .scale_cols(Wb), pa, pb),
               stats::cor(Wa[pa, ], Wb[pb, ]), tolerance = 1e-12)
})
