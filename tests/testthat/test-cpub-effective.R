test_that("CPUB-003 incompatible torch settings and frozen fixed mode fail early", {
  expect_error(.validate_backend_config("cebmf_torch", "point_laplace", 1,
    TRUE, TRUE, NULL, "residualize", list()), "distinct experimental")
  expect_error(.validate_backend_config("cebmf_torch", "laplace", 2,
    TRUE, FALSE, ~age, "fixed", list()), "no fixed-loadings")
  expect_error(.validate_backend_config("cebmf_torch", "laplace", 2,
    TRUE, FALSE, NULL, "residualize", list()), NA)
  expect_error(.validate_backend_config("flashier", "point_laplace", 1,
    TRUE, TRUE, NULL, "residualize", list(S = .1)), "explicit arguments")
})

test_that("CPUB-003 known noise changes actual ordinary fit noise and provenance", {
  set.seed(121)
  y <- matrix(rnorm(30 * 20), 30, 20,
              dimnames = list(paste0("d", 1:30), paste0("g", 1:20)))
  x <- as_cell_program_data(list(A = y))
  fit <- suppressWarnings(fit_celltype_programs(x, max_factors = 1,
    var_floor = .4, backfit = FALSE, nullcheck = FALSE))
  f <- fit$fits$A
  expect_equal(as.numeric(f$flash$residuals_sd), .4, tolerance = 1e-12)
  expect_identical(f$effective_config$noise_mode, "known_fixed")
  expect_null(f$effective_config$var_type)
  expect_equal(f$effective_config$noise_sd, .4)
  expect_identical(f$effective_config$loading_prior, "point_laplace")
  s <- freeze_preprocessing(x, covariates = NULL, feature_scale = "log1p_cp10k")
  f2 <- fit_celltype_programs(x, max_factors = 1, preprocessing_spec = s)
  expect_identical(f2$preprocessing$feature_scale, "log1p_cp10k")
  meta <- data.frame(observation_id = rownames(y), age = seq_len(30))
  x$sample_metadata <- meta
  s <- freeze_preprocessing(x, covariates = ~age, covariate_mode = "fixed")
  expect_error(fit_celltype_programs(x, max_factors = 1, preprocessing_spec = s,
                                    backend = "cebmf_torch"), "no fixed-loadings")
})
