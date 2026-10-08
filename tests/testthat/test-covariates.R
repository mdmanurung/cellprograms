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

test_that("a covariate constant within one cell type is skipped there, not an error", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  site <- ifelse(x$sample_metadata$batch == "a", "s1", "s2")
  x$sample_metadata$site <- site
  # cell type B only keeps donors from site s1 -> `site` is constant there
  keep <- rownames(x$pseudobulk$B)[site == "s1"]
  x$pseudobulk$B <- x$pseudobulk$B[keep, ]
  fit <- fit_celltype_programs(x, covariates = ~ site + age, max_factors = 2, seed = 1)
  XA <- cbind(1, site == "s2", x$sample_metadata$age)
  expect_lt(max(abs(crossprod(XA, fit$fits$A$Y_used))), 1e-8)     # A adjusted for site and age
  ageB <- x$sample_metadata$age[match(keep, x$sample_metadata$observation_id)]
  expect_lt(max(abs(crossprod(cbind(1, ageB), fit$fits$B$Y_used))), 1e-8)  # B: age only
})

test_that("covariate_mode = 'fixed' keeps covariate columns out of the programs and recovers the program", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  fit <- canonicalize_programs(fit_celltype_programs(x, covariates = ~ batch + age,
                                                     covariate_mode = "fixed", max_factors = 3, seed = 1))
  S <- fit$scores$A
  expect_gte(ncol(S), 1L)
  expect_identical(fit$preprocessing$covariate_mode, "fixed")
  fl <- fit$fits$A$flash
  expect_equal(ncol(fl$fixed_L), 2L)                       # batch (1 col) + age
  expect_lt(max(abs(crossprod(fl$fixed_L, S))) / nrow(S), 0.3)  # free programs ~ orthogonal to the covariates
  # Y_used is Y minus the fitted covariate part, so the summary still works
  expect_true(all(is.finite(program_summary(fit)$R2)))
  expect_error(fit_celltype_programs(x, covariates = ~ age, covariate_mode = "fixed", center = FALSE), "center")
})
