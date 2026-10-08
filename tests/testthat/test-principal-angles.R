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
