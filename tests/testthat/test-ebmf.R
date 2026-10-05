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
  fit <- canonicalize_programs(fit_celltype_programs(cpd, seed = 1))
  expect_true(isTRUE(fit$fits$B$skipped))
  expect_equal(ncol(fit$scores$B), 0)
  expect_equal(nrow(fit$scores$B), 10)
  expect_equal(ncol(fit$scores$A) >= 0, TRUE)
})
