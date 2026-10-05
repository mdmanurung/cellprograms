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
