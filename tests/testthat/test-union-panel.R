# C5 route-B helper: post-fit union-panel loadings on fixed scores.

.fake_fit <- function(Y, S, preprocessing, meta, fitted_genes = colnames(Y)) {
  fake <- list(
    fits = list(A = list(genes = fitted_genes)),
    scores = list(A = S),
    pseudobulk = list(A = Y),
    preprocessing = preprocessing,
    sample_metadata = meta
  )
  class(fake) <- c("cell_program_fit", "cell_program_data", "list")
  fake
}

test_that("union_panel_loadings: least-squares refit matches hand-computed replay (residualize)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(21)
  n <- 12; g <- 4
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", seq_len(n)), paste0("g", seq_len(g))))
  batch <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
  Y[, 1] <- 3 * batch                     # pure batch gene: must get ~0 loading after replay
  S <- matrix(c(-1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0), ncol = 1,
              dimnames = list(rownames(Y), "A_1"))
  meta <- data.frame(observation_id = rownames(Y), batch = factor(batch))
  pp <- list(centered = TRUE, scaled = FALSE, covariates = ~batch, covariate_mode = "residualize")
  fake <- .fake_fit(Y, S, pp, meta, fitted_genes = paste0("g", seq_len(g)))

  up <- union_panel_loadings(fake)
  expect_setequal(up$genes, paste0("g", seq_len(g)))
  X <- cbind(1, batch)
  Yres <- stats::lm.fit(X, Y)$residuals
  expected <- t(stats::lm.fit(S, Yres)$coefficients)
  expect_equal(up$loadings$A[paste0("g", seq_len(g)), , drop = FALSE], expected)
  expect_equal(unname(up$loadings$A["g1", 1]), 0, tolerance = 1e-12)
  # discrimination: without the residualization replay the batch gene loads strongly
  no_replay <- t(stats::lm.fit(S, scale(Y))$coefficients)
  expect_true(abs(no_replay["g1", 1]) > 0.3)
})

test_that("union_panel_loadings: fixed-mode replay subtracts the OLS covariate part", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(22)
  n <- 12; g <- 4
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", seq_len(n)), paste0("g", seq_len(g))))
  batch <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
  Y[, 1] <- 3 * batch
  S <- matrix(c(-1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0), ncol = 1,
              dimnames = list(rownames(Y), "A_1"))
  meta <- data.frame(observation_id = rownames(Y), batch = factor(batch))
  Xc <- matrix(batch - mean(batch), ncol = 1, dimnames = list(rownames(Y), "batchb"))
  pp <- list(centered = TRUE, scaled = FALSE, covariates = ~batch, covariate_mode = "fixed")
  fake <- .fake_fit(Y, S, pp, meta, fitted_genes = paste0("g", seq_len(g)))
  fake$fits$A$flash <- list(fixed_L = Xc)

  up <- union_panel_loadings(fake)
  F0 <- solve(crossprod(Xc), crossprod(Xc, Y))
  expected <- t(stats::lm.fit(S, scale(Y - Xc %*% F0, center = TRUE, scale = FALSE))$coefficients)
  expect_equal(up$loadings$A[paste0("g", seq_len(g)), , drop = FALSE], expected)
  expect_equal(unname(up$loadings$A["g1", 1]), 0, tolerance = 1e-12)
})

test_that("union_panel_loadings: end-to-end on a features='variable' fit", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, features = "variable",
                                                     n_variable_genes = 150, seed = 1))
  up <- union_panel_loadings(fit)
  nr <- sapply(up$loadings, nrow)
  expect_true(all(nr == length(up$genes)))
  expect_true(length(up$genes) > 150)   # union exceeds any single type's panel
  expect_identical(rownames(up$loadings[[1]]), rownames(up$loadings[[2]]))
  for (ct in names(fit$fits)) {
    fitted <- fit$fits[[ct]]$genes
    expect_true(all(is.finite(up$loadings[[ct]][fitted, ])))
    # least squares on fixed scores is optimal: refit R2 >= flash-loading R2
    S <- fit$scores[[ct]]
    Y <- fit$fits[[ct]]$Y_used
    res <- Y - S %*% t(up$loadings[[ct]][fitted, , drop = FALSE])
    r2_refit <- 1 - sum(res^2) / sum(sweep(Y, 2, colMeans(Y))^2)
    expect_gte(r2_refit, program_summary(fit)$R2[match(ct, program_summary(fit)$cell_type)] - 1e-6)
    # refit loadings track the shrunk flashier posterior means on fitted genes
    for (k in seq_len(ncol(S))) {
      cr <- suppressWarnings(cor(up$loadings[[ct]][fitted, k], fit$loadings[[ct]][fitted, k]))
      if (is.finite(cr)) expect_gt(cr, 0.9)
    }
  }
})

test_that("union_panel_loadings: skipped cell types, custom panels, and input validation", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(7)
  g <- 20
  Y1 <- matrix(rnorm(10 * g), 10, g, dimnames = list(paste0("s", seq_len(10)), paste0("g", seq_len(g))))
  Y1[, 1:5] <- Y1[, 1:5] + outer(rnorm(10), rep(1, 5)) * 2   # plant one program: K >= 1
  Y2 <- matrix(rep(rnorm(g), each = 10), 10, g, dimnames = list(paste0("s", seq_len(10)), paste0("g", seq_len(g))))
  cpd <- as_cell_program_data(list(A = Y1, B = Y2))
  fit <- suppressMessages(canonicalize_programs(
    fit_celltype_programs(cpd, max_factors = 3, seed = 1)))
  up <- union_panel_loadings(fit)
  expect_equal(ncol(up$loadings$B), 0L)          # degenerate type: K = 0
  expect_equal(nrow(up$loadings$B), length(up$genes))
  expect_gt(ncol(up$loadings$A), 0L)

  up2 <- union_panel_loadings(fit, genes = c("g1", "not_a_gene"))
  expect_identical(up2$genes, c("g1", "not_a_gene"))
  expect_equal(unname(up2$loadings$A["not_a_gene", 1]), NA_real_)

  expect_error(union_panel_loadings(list(a = 1)), "cell_program_fit")
  expect_error(union_panel_loadings(fit, genes = character(0)), "non-empty")
  expect_error(union_panel_loadings(fit, genes = 3), "character")
})
