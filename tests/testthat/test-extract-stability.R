test_that("extract helpers return labelled, consistent objects", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  fit <- canonicalize_programs(fit_celltype_programs(x, covariates = ~ batch, max_factors = 3, seed = 1))
  meta <- program_metadata(fit)
  Z <- program_scores(fit)
  expect_equal(colnames(Z), meta$program_id)
  expect_equal(rownames(Z), x$observation_ids)
  expect_equal(nrow(program_scores(fit, "long")), sum(vapply(fit$scores, length, 1L)))
  expect_equal(program_summary(fit)$n_programs, as.integer(table(meta$cell_type)[c("A", "B")]))
  tg <- top_genes(fit, meta$program_id[1], n = 5, by = "lfsr")
  expect_equal(nrow(tg), 5L)
})

test_that("assess_program_stability recovers the strong planted program", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), covariates = ~ batch + age,
                                                     max_factors = 3, seed = 1))
  fit <- assess_program_stability(fit, n_boot = 6, seed = 1, progress = FALSE)
  st <- fit$stability$per_program
  expect_equal(st$program_id, program_metadata(fit)$program_id)
  expect_gt(max(st$recovery_freq), 0.8)
})

test_that("stability refits keep covariate_mode and flash_control", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), covariates = ~ batch,
                                                     covariate_mode = "fixed", max_factors = 3,
                                                     flash_control = list(backfit = FALSE), seed = 1))
  seen <- list()
  orig <- fit_celltype_programs
  ## mock inside the package namespace; assigning into globalenv is invisible to package code
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...); seen[[length(seen) + 1L]] <<- a[c("covariate_mode", "flash_control")]
    orig(...)
  })
  assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE)
  expect_length(seen, 2L)
  for (s in seen) {
    expect_identical(s$covariate_mode, "fixed")
    expect_identical(s$flash_control, list(backfit = FALSE))
  }
})

test_that("fit consumers fail closed with a clear message on a non-canonicalized fit", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 200, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  fit_raw <- fit_celltype_programs(x, max_factors = 2, seed = 1)  # NOT canonicalized
  expect_null(fit_raw$canonicalization)

  guarded <- list(
    ~program_summary(fit_raw),
    ~program_metadata(fit_raw),
    ~program_scores(fit_raw),
    ~program_loadings(fit_raw),
    ~top_genes(fit_raw, "B_1"),
    ~assess_program_stability(fit_raw, n_boot = 2, seed = 1),
    ~principal_angles(fit_raw, "B", "CD4", n_perm = 20),
    ~sharing_spectrum(fit_raw, n_perm = 20),
    ~align_programs(fit_raw),
    ~match_programs(fit_raw, fit_raw, "B")
  )
  for (g in guarded) {
    expect_error(eval(g[[2]], envir = environment(g)), "canonicalize_programs",
                 label = "raw fit rejected with the fix named")
  }
  # and the guard is silent on a canonicalized fit
  fit_c <- canonicalize_programs(fit_raw)
  expect_s3_class(program_summary(fit_c), "data.frame")
})
