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

test_that("a flashier nullcheck failure falls back to no nullcheck and is recorded", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  real <- flashier::flash
  testthat::local_mocked_bindings(
    flash = function(...) { if (isTRUE(list(...)$nullcheck)) stop("NaN ELBO") else real(...) },
    .package = "flashier")
  fit <- suppressWarnings(fit_celltype_programs(.toy(), max_factors = 2,
                                                loading_prior = "point_normal", seed = 1))
  expect_identical(fit$fits$A$fallback, "no_nullcheck")
  expect_warning(fit_celltype_programs(.toy(), max_factors = 2, loading_prior = "point_normal"), "fallback")
})

## Regression tests for the 2026-10 covariate audit (reports/audit/covariate_audit/):
## the four confirmed defects, asserted at their FIXED behavior.

test_that("F1: replay center/scale are subset when zero-variance genes drop (default path)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  x$pseudobulk$A[, "g1"] <- 5  # constant in A -> zero-variance after centering
  msgs <- testthat::capture_messages(
    fit <- fit_celltype_programs(x, max_factors = 2, features = "all",
                                 scale = TRUE, feature_scale = "logcounts", seed = 1))
  expect_true(any(grepl("zero-variance", msgs)))  # the drop is announced
  rp <- fit$fits$A$preproc_replay
  genes <- fit$fits$A$genes
  expect_false("g1" %in% genes)
  expect_length(rp$center, length(genes))   # was 60 (full panel) before the fix
  expect_length(rp$scale, length(genes))
  # aligned replay reproduces Y_used exactly (misaligned by >10 before the fix)
  Q <- as.matrix(x$pseudobulk$A)[, genes, drop = FALSE]
  Qc <- Q - matrix(rp$center, nrow = nrow(Q), ncol = length(genes), byrow = TRUE)
  if (!is.null(rp$scale)) Qc <- Qc / matrix(rp$scale, nrow = nrow(Q), ncol = length(genes), byrow = TRUE)
  expect_lt(max(abs(Qc - fit$fits$A$Y_used[, genes])), 1e-8)
})

test_that("F2: collinear design fits and projects (NA aliased coefficients treated as 0)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  x$sample_metadata$age_copy <- x$sample_metadata$age
  fit <- fit_celltype_programs(x, covariates = c("age", "age_copy"), max_factors = 2,
                               feature_scale = "logcounts", seed = 1)
  beta <- fit$fits$A$preproc_replay$residualize_coefficients
  expect_true(anyNA(beta))                       # lm.fit pivot: aliased column is NA
  fit_k <- canonicalize_programs(fit)
  ids_a <- rownames(x$pseudobulk$A)
  pr <- project_programs_to_cells(
    fit_k, x$pseudobulk$A, cell_types = setNames(rep("A", length(ids_a)), ids_a),
    donors = setNames(ids_a, ids_a), sample_metadata = x$sample_metadata)
  expect_true(all(is.finite(pr$activity[, 1])))
})

test_that("F3: a frozen spec with covariate_mode='fixed', center=FALSE is refused", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  spec_bad <- freeze_preprocessing(x, features = "all", covariates = ~ age,
                                   covariate_mode = "fixed", center = FALSE)
  expect_error(fit_celltype_programs(x, max_factors = 2, preprocessing_spec = spec_bad),
               "center")
})

test_that("F4: interaction referencing a within-ct-constant covariate is dropped cleanly", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  x$sample_metadata$sex <- factor(rep(c("F", "M"), each = 20))
  x$pseudobulk$A <- x$pseudobulk$A[1:20, ]   # sex constant (all F) within A
  expect_message(
    fit <- fit_celltype_programs(x, covariates = ~ age + sex + age:sex,
                                 max_factors = 2, seed = 1),
    "constant in cell\\s+type")
  # A adjusted for age only (sex + interaction dropped there); still orthogonal to age
  ageA <- x$sample_metadata$age[1:20]
  expect_lt(max(abs(crossprod(cbind(1, ageA), fit$fits$A$Y_used))), 1e-8)
})
