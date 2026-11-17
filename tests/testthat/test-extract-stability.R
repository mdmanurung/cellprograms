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
  # Whether a refit needs the point_normal fallback depends on the installed
  # flashier/ebnm, so the test checks that the log records one fallback entry per
  # cell type and refit, not that the fallback fires.
  fit <- suppressWarnings(assess_program_stability(fit, n_boot = 6, seed = 1, progress = FALSE))
  effective <- fit$stability$replicate_log$effective_refit
  recorded <- unlist(lapply(effective, function(refit) {
    vapply(refit, function(ct) as.character(ct$fallback)[1], character(1))
  }), use.names = FALSE)
  expect_length(recorded, sum(lengths(effective)))
  expect_false(anyNA(recorded))
  expect_identical(fit$stability$target,
                   "full_fitting_procedure_including_documented_fallback_ladder")
  st <- fit$stability$per_program
  expect_equal(st$program_id, program_metadata(fit)$program_id)
  expect_gt(max(st$recovery_freq), 0.8)
})

test_that("stability refits keep covariate_mode and explicit backfit control", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), covariates = ~ batch,
                                                     covariate_mode = "fixed", max_factors = 3,
                                                     backfit = FALSE, seed = 1))
  seen <- list()
  orig <- fit_celltype_programs
  ## mock inside the package namespace; assigning into globalenv is invisible to package code
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...); seen[[length(seen) + 1L]] <<- a[c("covariate_mode", "backfit")]
    orig(...)
  })
  assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE)
  expect_length(seen, 2L)
  for (s in seen) {
    expect_identical(s$covariate_mode, "fixed")
    expect_identical(s$backfit, FALSE)
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

test_that("match_method is deprecated and does not change stability results", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), covariates = ~ batch + age,
                                                     max_factors = 3, seed = 1))
  run <- function(...) suppressWarnings(assess_program_stability(
    fit, n_boot = 3, seed = 1, progress = FALSE, ...))
  base <- run()
  expect_warning(sub <- assess_program_stability(fit, n_boot = 3, seed = 1, progress = FALSE,
                                                 match_method = "subspace"),
                 "match_method.*deprecated")
  corr <- run(match_method = "correlation")
  expect_equal(sub$stability$per_program, base$stability$per_program)
  expect_equal(corr$stability$per_program, base$stability$per_program)
})

# ---- repeated visits: donor-level draws ---------------------------------------------------

.toy_visits <- function() {
  x <- .toy(n = 40)
  x$sample_metadata$donor <- rep(paste0("d", 1:20), each = 2)
  x$sample_metadata$visit <- rep(c("v1", "v2"), 20)
  x$sample_metadata$site <- rep(c("s1", "s2"), each = 20)
  x
}

test_that("default stability draws keep the RNG stream of the pre-change code at seed = 1", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- .pa_golden_fit()
  draws <- list()
  orig <- fit_celltype_programs
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...)
    draws[[length(draws) + 1L]] <<- lapply(a[[1]]$pseudobulk, \(Y) match(rownames(Y), paste0("obs", 1:40)))
    orig(...)
  })
  suppressWarnings(assess_program_stability(fit, n_boot = 3, seed = 1, progress = FALSE))
  # Sampled row positions per refit and cell type, captured from release/0.3.0 @ 166e40c
  # (before the donor-draw change). They do not depend on flashier or ebnm numerics.
  expected <- list(list(B = c(4L, 39L, 1L, 34L, 23L, 14L, 18L, 33L, 21L, 32L, 
  10L, 22L, 35L, 30L, 7L, 9L, 15L, 31L, 5L, 25L, 28L, 29L, 19L, 
  2L, 27L, 20L, 36L, 12L, 38L, 40L, 3L, 6L), CD4 = c(10L, 38L, 
  20L, 28L, 39L, 23L, 6L, 25L, 35L, 24L, 34L, 26L, 7L, 19L, 40L, 
  30L, 31L, 14L, 2L, 13L, 18L, 32L, 33L, 1L, 3L, 16L, 8L, 36L, 
  37L, 21L, 27L, 12L), CD8 = c(11L, 17L, 36L, 13L, 25L, 38L, 23L, 
  20L, 7L, 29L, 37L, 28L, 22L, 12L, 16L, 1L, 30L, 21L, 6L, 39L, 
  9L, 32L, 18L, 26L, 40L, 14L, 10L, 19L, 27L, 3L, 2L, 35L), NK = c(26L, 
  15L, 29L, 24L, 12L, 37L, 28L, 22L, 8L, 1L, 19L, 40L, 3L, 11L, 
  31L, 14L, 33L, 6L, 20L, 9L, 5L, 25L, 17L, 34L, 23L, 18L, 7L, 
  39L, 13L, 38L, 21L, 2L), Mono = c(27L, 7L, 20L, 24L, 16L, 19L, 
  2L, 1L, 12L, 39L, 9L, 28L, 26L, 13L, 29L, 34L, 25L, 35L, 17L, 
  37L, 11L, 3L, 23L, 14L, 18L, 40L, 38L, 8L, 30L, 21L, 10L, 15L
  )), list(B = c(10L, 15L, 37L, 32L, 31L, 29L, 12L, 14L, 28L, 6L, 
  20L, 16L, 39L, 2L, 1L, 3L, 24L, 9L, 36L, 17L, 30L, 5L, 23L, 34L, 
  38L, 19L, 8L, 4L, 26L, 18L, 35L, 21L), CD4 = c(20L, 24L, 22L, 
  16L, 30L, 29L, 2L, 8L, 36L, 14L, 4L, 26L, 33L, 21L, 5L, 28L, 
  39L, 7L, 40L, 32L, 17L, 10L, 27L, 9L, 11L, 38L, 15L, 6L, 13L, 
  18L, 3L, 19L), CD8 = c(13L, 7L, 25L, 31L, 28L, 36L, 20L, 29L, 
  37L, 10L, 35L, 8L, 1L, 14L, 12L, 16L, 26L, 32L, 2L, 17L, 30L, 
  4L, 5L, 24L, 15L, 27L, 34L, 22L, 38L, 33L, 18L, 11L), NK = c(39L, 
  7L, 27L, 9L, 20L, 34L, 29L, 11L, 2L, 36L, 38L, 3L, 15L, 6L, 37L, 
  28L, 23L, 5L, 14L, 24L, 33L, 31L, 17L, 1L, 10L, 13L, 18L, 35L, 
  16L, 40L, 4L, 22L), Mono = c(27L, 39L, 38L, 16L, 34L, 29L, 10L, 
  6L, 28L, 7L, 17L, 1L, 23L, 35L, 3L, 11L, 13L, 32L, 18L, 5L, 37L, 
  8L, 2L, 4L, 19L, 9L, 33L, 40L, 20L, 15L, 25L, 14L)), list(B = c(32L, 
  3L, 19L, 33L, 4L, 13L, 1L, 5L, 10L, 12L, 39L, 25L, 23L, 21L, 
  30L, 36L, 9L, 28L, 26L, 22L, 20L, 24L, 31L, 18L, 7L, 35L, 11L, 
  29L, 27L, 16L, 34L, 37L), CD4 = c(25L, 9L, 21L, 38L, 15L, 24L, 
  14L, 7L, 32L, 27L, 31L, 26L, 13L, 29L, 36L, 18L, 6L, 2L, 37L, 
  1L, 17L, 28L, 19L, 8L, 23L, 39L, 30L, 20L, 40L, 34L, 4L, 10L), 
      CD8 = c(5L, 10L, 20L, 15L, 29L, 2L, 1L, 40L, 26L, 25L, 33L, 
      8L, 13L, 31L, 4L, 39L, 6L, 18L, 30L, 38L, 34L, 24L, 14L, 
      28L, 36L, 21L, 3L, 12L, 37L, 23L, 19L, 16L), NK = c(23L, 
      20L, 5L, 2L, 38L, 12L, 3L, 26L, 18L, 14L, 13L, 40L, 15L, 
      24L, 37L, 30L, 16L, 7L, 27L, 17L, 11L, 9L, 8L, 25L, 6L, 34L, 
      39L, 35L, 10L, 1L, 22L, 31L), Mono = c(13L, 26L, 2L, 32L, 
      36L, 5L, 23L, 20L, 28L, 37L, 12L, 3L, 4L, 15L, 14L, 40L, 
      29L, 39L, 27L, 38L, 6L, 19L, 16L, 31L, 18L, 10L, 35L, 34L, 
      7L, 22L, 25L, 8L)))
  expect_identical(draws, expected)
  # seed = 1L is the default
  draws <- list()
  suppressWarnings(assess_program_stability(fit, n_boot = 3, progress = FALSE))
  expect_identical(draws, expected)
})

test_that("donor_col draws whole donors, the same donors for every cell type", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy_visits()
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 2, seed = 1))
  draws <- list()
  orig <- fit_celltype_programs
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...)
    draws[[length(draws) + 1L]] <<- lapply(a[[1]]$pseudobulk, rownames)
    orig(...)
  })
  st <- suppressWarnings(assess_program_stability(fit, n_boot = 3, seed = 2, progress = FALSE,
                                                  donor_col = "donor", visit_col = "visit"))
  expect_length(draws, 3L)
  donor_of <- stats::setNames(x$sample_metadata$donor, x$sample_metadata$observation_id)
  for (d in draws) {
    donors_by_ct <- lapply(d, function(ids) sort(unique(donor_of[ids])))
    expect_identical(donors_by_ct$A, donors_by_ct$B)                    # one donor set for all cell types
    expect_length(donors_by_ct$A, 16L)                                  # round(0.8 * 20)
    expect_true(all(table(donor_of[d$A]) == 2L))                        # no donor split across visits
  }
  # same seed, same draws
  draws2 <- list()
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...); draws2[[length(draws2) + 1L]] <<- lapply(a[[1]]$pseudobulk, rownames); orig(...)
  })
  suppressWarnings(assess_program_stability(fit, n_boot = 3, seed = 2, progress = FALSE,
                                            donor_col = "donor", visit_col = "visit"))
  expect_identical(draws2, draws)
})

test_that("donor_col strata draw proportionally and the design is validated", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy_visits()
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 2, seed = 1))
  expect_error(assess_program_stability(fit, n_boot = 1, progress = FALSE, donor_col = "donor"),
               "visit_col")
  expect_error(assess_program_stability(fit, n_boot = 1, progress = FALSE, donor_col = "donor",
                                        visit_col = "visit", donor_strata = "visit"),
               "constant within donor")
  expect_error(assess_program_stability(fit, n_boot = 1, progress = FALSE, visit_col = "visit"),
               "requires `donor_col`")
  draws <- list()
  orig <- fit_celltype_programs
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...); draws[[length(draws) + 1L]] <<- rownames(a[[1]]$pseudobulk$A); orig(...)
  })
  suppressWarnings(assess_program_stability(fit, n_boot = 2, seed = 3, progress = FALSE,
                                            donor_col = "donor", visit_col = "visit",
                                            donor_strata = "site"))
  site_of <- stats::setNames(x$sample_metadata$site, x$sample_metadata$observation_id)
  for (ids in draws) expect_equal(as.integer(table(site_of[ids])), c(16L, 16L))   # 8 donors x 2 visits per site
})
