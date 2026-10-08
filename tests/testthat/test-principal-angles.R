test_that("principal_angles on adjusted fits finds the shared program", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), covariates = c("batch", "age"),
                                                     max_factors = 3, seed = 1))
  pa <- principal_angles(fit, "A", "B", n_perm = 99, seed = 1)
  expect_gt(pa$cosines[1], 0.9)
  expect_lt(pa$p_values[1], 0.05)
  pl <- principal_angles(fit, "A", "B", space = "loadings", weight = "residual_sd", n_perm = 0)
  expect_s3_class(pl, "principal_angles")
})

test_that(".perm_within permutes only inside strata", {
  set.seed(1)
  st <- rep(c("x", "y", "z"), c(5, 6, 7))
  for (i in 1:20) {
    p <- .perm_within(length(st), st)
    expect_identical(st[p], st)
    expect_setequal(p, seq_along(st))
  }
})

test_that("strata is validated and accepted", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  ids <- intersect(rownames(fit$scores$A), rownames(fit$scores$B))
  st <- setNames(rep(c("u", "v"), length.out = length(ids)), ids)
  pa <- principal_angles(fit, "A", "B", n_perm = 19, seed = 1, strata = st)
  expect_s3_class(pa, "principal_angles")
  expect_error(principal_angles(fit, "A", "B", n_perm = 5, strata = unname(st)), "named vector")
  expect_error(principal_angles(fit, "A", "B", space = "loadings", n_perm = 5, strata = st), "scores space")
})
