test_that("covariates= residualizes Y_used orthogonally to the covariates", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  fit <- fit_celltype_programs(x, covariates = ~ batch + age, max_factors = 3, seed = 1)
  Yu <- fit$fits$A$Y_used
  X <- cbind(1, x$sample_metadata$batch == "b", x$sample_metadata$age)
  expect_lt(max(abs(crossprod(X, Yu))), 1e-8)
  expect_identical(fit$preprocessing$covariates, ~ batch + age)
  expect_error(fit_celltype_programs(x, covariates = "nope"), "not in")
})
