test_that("fit_celltype_programs runs flashier end-to-end and canonicalizes", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  set.seed(11)
  n <- 30; g <- 60
  Y <- matrix(rnorm(n * g), n, g)
  Y[1:15, 1:5] <- Y[1:15, 1:5] + 3  # planted program in one sample group
  rownames(Y) <- paste0("s", seq_len(n))
  colnames(Y) <- paste0("g", seq_len(g))
  cpd <- as_cell_program_data(list(CT1 = Y),
                              sample_metadata = data.frame(observation_id = rownames(Y)))
  fit <- canonicalize_programs(fit_celltype_programs(cpd, seed = 1))
  S <- fit$scores$CT1
  W <- fit$loadings$CT1
  expect_true(is.matrix(S))
  expect_equal(rownames(S), rownames(Y))
  expect_equal(nrow(W), g)
  expect_equal(ncol(S), ncol(W))
  # canonicalization preserves the raw flashier reconstruction exactly
  raw <- fit$fits$CT1$flash
  expect_lt(max(abs(tcrossprod(S, W) -
                     tcrossprod(as.matrix(raw$L_pm), as.matrix(raw$F_pm)))), 1e-8)
})

test_that("fully degenerate cell type yields K=0 fit instead of error", {
  set.seed(7)
  g <- 20
  Y1 <- matrix(rnorm(10 * g), 10, g, dimnames = list(paste0("s", 1:10), paste0("g", 1:g)))
  # cell type B: every row identical -> all genes constant -> zero after centering
  Y2 <- matrix(rep(rnorm(g), each = 10), 10, g, dimnames = list(paste0("s", 1:10), paste0("g", 1:g)))
  cpd <- as_cell_program_data(list(A = Y1, B = Y2))
  # max_factors = 3: the C4 residual-df guard caps K at n - q - 1 = 9 here, and
  # this test exercises the degenerate sentinel, not the K:n ratio.
  fit <- canonicalize_programs(fit_celltype_programs(cpd, max_factors = 3, seed = 1))
  expect_true(isTRUE(fit$fits$B$skipped))
  expect_equal(ncol(fit$scores$B), 0)
  expect_equal(nrow(fit$scores$B), 10)
  expect_equal(ncol(fit$scores$A) >= 0, TRUE)
})


test_that("var_type validation rejects invalid specs", {
  skip_if_not_installed("flashier")
  g <- 20
  Y <- matrix(rnorm(10 * g), 10, g, dimnames = list(paste0("s", 1:10), paste0("g", 1:g)))
  cpd <- as_cell_program_data(list(CT1 = Y))
  expect_error(fit_celltype_programs(cpd, var_type = 3), "var_type")
  expect_error(fit_celltype_programs(cpd, var_type = c(0, 1)), "var_type")
  expect_error(fit_celltype_programs(cpd, var_type = "1"), "var_type")
})

test_that("var_type = c(1, 2) Kronecker fits end-to-end and is recorded", {
  skip_if_not_installed("flashier")
  n <- 12
  g <- 30
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", seq_len(n)), paste0("g", seq_len(g))))
  Y[1:6, 1:5] <- Y[1:6, 1:5] + 3
  cpd <- as_cell_program_data(list(CT1 = Y),
                              sample_metadata = data.frame(observation_id = rownames(Y)))
  fit <- canonicalize_programs(fit_celltype_programs(cpd, var_type = c(1, 2), seed = 1))
  S <- fit$scores$CT1
  W <- fit$loadings$CT1
  expect_true(is.matrix(S))
  expect_equal(rownames(S), rownames(Y))
  expect_equal(ncol(S), ncol(W))
  expect_equal(fit$preprocessing$var_type, c(1, 2))
})

test_that("C4 guard: max_factors above the residual degrees of freedom aborts", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(11)
  g <- 30
  mk <- function(n, covars = NULL) {
    Y <- matrix(rnorm(n * g, sd = 0.5), n, g,
                dimnames = list(paste0("s", seq_len(n)), paste0("g", seq_len(g))))
    meta <- data.frame(observation_id = rownames(Y))
    if (!is.null(covars)) meta <- cbind(meta, covars)
    as_cell_program_data(list(A = Y), sample_metadata = meta)
  }
  # No covariates (where no other guard protects): n - q - 1 = 7 < 10 factors.
  expect_error(fit_celltype_programs(mk(8), max_factors = 10, seed = 1),
               "residual degrees of freedom")
  # Covariate columns count toward q: n = 14, q = 2 -> bound 11 < 12 factors,
  # while the older 10-residual-df guard (n - rk = 11) still passes.
  meta <- data.frame(batch = factor(rep(c("a", "b"), length.out = 14)),
                     age = rnorm(14))
  expect_error(fit_celltype_programs(mk(14, meta), covariates = ~ batch + age,
                                     max_factors = 12, seed = 1),
               "residual degrees of freedom")
  # Fixed mode counts the fixed columns: q = 1 -> bound 12 < 13 factors.
  expect_error(fit_celltype_programs(mk(14, meta), covariates = ~batch,
                                     covariate_mode = "fixed", max_factors = 13, seed = 1),
               "residual degrees of freedom")
})

test_that("C4 guard does not bind at or below the bound, and degenerate types keep precedence", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(12)
  g <- 30
  mk <- function(Y_list, meta = NULL) {
    as_cell_program_data(Y_list, sample_metadata = meta)
  }
  # Exactly at the bound: n - q - 1 = 11 >= 11 factors -> fits.
  Y <- matrix(rnorm(12 * g, sd = 0.5), 12, g,
              dimnames = list(paste0("s", seq_len(12)), paste0("g", seq_len(g))))
  cpd <- mk(list(A = Y), data.frame(observation_id = rownames(Y)))
  fit <- fit_celltype_programs(cpd, max_factors = 11, loading_prior = "point_normal", seed = 1)
  expect_false(isTRUE(fit$fits$A$skipped))
  # Well below the bound on few donors: n = 8, K = 2 fits.
  Y8 <- matrix(rnorm(8 * g, sd = 0.5), 8, g,
               dimnames = list(paste0("s", seq_len(8)), paste0("g", seq_len(g))))
  cpd8 <- mk(list(A = Y8), data.frame(observation_id = rownames(Y8)))
  fit8 <- fit_celltype_programs(cpd8, max_factors = 2, loading_prior = "point_normal", seed = 1)
  expect_false(isTRUE(fit8$fits$A$skipped))
  # A fully degenerate cell type records K = 0 even when max_factors would violate C4
  # (the degenerate sentinel is checked before the guard, and no real fit is attempted).
  Yc <- matrix(rep(rnorm(g), each = 8), 8, g,
               dimnames = list(paste0("s", seq_len(8)), paste0("g", seq_len(g))))
  cpdd <- mk(list(A = Yc), data.frame(observation_id = rownames(Yc)))
  fitd <- suppressMessages(fit_celltype_programs(cpdd, max_factors = 10,
                                                 loading_prior = "point_normal", seed = 1))
  expect_true(isTRUE(fitd$fits$A$skipped))
})
