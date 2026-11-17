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
  fit <- .expect_warning_messages(
    assess_program_stability(fit, n_boot = 6, seed = 1, progress = FALSE),
    'flashier needed fallback "point_normal" in cell type "B".')
  effective <- fit$stability$replicate_log$effective_refit
  recorded <- unlist(lapply(effective, function(refit) {
    vapply(refit, function(ct) ct$fallback, character(1))
  }), use.names = FALSE)
  expect_true(any(recorded == "point_normal"))
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

test_that("default stability results keep the RNG stream of the pre-change code at seed = 1", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- .pa_golden_fit()
  # golden values captured from release/0.3.0 @ 166e40c, before this change
  st <- suppressWarnings(assess_program_stability(fit, n_boot = 3, seed = 1, progress = FALSE))$stability$per_program
  expect_equal(st$median_score_cor,
               c(0.9999290303, 0.9996772587, 0.9999642101, 0.999873038, 0.9999572966,
                 0.9999132847, 0.9999318474, 0.9999165394), tolerance = 1e-8)
  expect_equal(st$top_gene_jaccard,
               c(1, 1, 0.9047619048, 0.8181818182, 0.9047619048, 0.9047619048, 1, 0.9047619048),
               tolerance = 1e-8)
  # seed = 1L is now the default
  st0 <- suppressWarnings(assess_program_stability(fit, n_boot = 3, progress = FALSE))$stability$per_program
  expect_equal(st0, st)
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
