# Fixes from the external scientific-correctness review (Biomni, 2026-10-08)

test_that("fixed-mode design survives a factor with an unused level", {
  meta <- data.frame(observation_id = paste0("s", 1:20),
                     site = factor(rep(c("a", "b"), 10), levels = c("a", "b", "c")))
  X <- .cov_design(meta$observation_id, meta, "site", "CT", intercept = FALSE)
  expect_equal(ncol(X), 1L)
  expect_false(anyNA(X))
})

test_that("a level with a single donor triggers a leverage warning; too few residual df aborts", {
  meta <- data.frame(observation_id = paste0("s", 1:30),
                     g = c(rep("a", 14), rep("b", 15), "c"))
  expect_warning(.cov_design(meta$observation_id, meta, "g", "CT"), "leverage")
  small <- data.frame(observation_id = paste0("s", 1:12), g = rep(letters[1:6], 2))
  expect_error(.cov_design(small$observation_id, small, "g", "CT"), "residual degrees of freedom")
})

test_that("per-pair seeds differ between pairs and are reproducible", {
  expect_false(.pair_seed(1, c("A", "B")) == .pair_seed(1, c("A", "C")))
  expect_false(.pair_seed(1, c("A", "B")) == .pair_seed(1, c("B", "A")))
  expect_identical(.pair_seed(1, c("A", "B")), .pair_seed(1, c("A", "B")))
})

test_that("sharing_spectrum scales n_perm with the number of pairs and warns when too small", {
  set.seed(2); n <- 30
  mk <- function(ct) matrix(rnorm(n * 2), n, 2, dimnames = list(paste0("s", 1:n), paste0(ct, "_", 1:2)))
  sc <- list(A = mk("A"), B = mk("B"), C = mk("C"))
  fit <- structure(list(scores = sc, loadings = lapply(sc, function(S) S[0, ])), class = "cell_program_fit")
  expect_no_warning(r <- sharing_spectrum(fit, seed = 1))
  expect_gte(r$pairs[[1]]$n_perm, 999L)
  expect_warning(sharing_spectrum(fit, n_perm = 50, seed = 1), "permutation floor")
})

test_that("gene weights are capped at 10x the median", {
  set.seed(3)
  Y <- matrix(rnorm(40 * 6), 40, 6, dimnames = list(paste0("s", 1:40), paste0("g", 1:6)))
  Z <- matrix(rnorm(40), 40, 1, dimnames = list(rownames(Y), "A_1"))
  W <- matrix(c(1, rep(0, 5)), 6, 1, dimnames = list(colnames(Y), "A_1"))
  Y[, 1] <- Z[, 1]                         # gene 1 is fitted perfectly -> residual sd ~ 0
  fit <- list(fits = list(A = list(Y_used = Y), B = list(Y_used = Y)),
              scores = list(A = Z, B = Z), loadings = list(A = W, B = W))
  w <- .gene_weights(fit, "A", "B", colnames(Y))
  expect_lte(max(w) / median(w), 10 + 1e-8)
})

test_that("a silent NaN ELBO triggers the fallback ladder", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  real <- flashier::flash
  testthat::local_mocked_bindings(
    flash = function(...) { f <- real(...); if (isTRUE(list(...)$nullcheck)) f$elbo <- NaN; f },
    .package = "flashier")
  fit <- suppressWarnings(fit_celltype_programs(.toy(), max_factors = 2,
                                                loading_prior = "point_normal", seed = 1))
  expect_identical(fit$fits$A$fallback, "no_nullcheck")
})

test_that("degenerate permutation null gives NA z instead of NaN/Inf", {
  one <- matrix(1, 20, 1)   # constant column: every row permutation gives the same statistic, sd(null) = 0
  r <- .pa_engine(one, one, n_perm = 20, seed = 1)
  expect_true(is.na(r$z_pair))
  expect_true(is.finite(r$p_pair))
})
