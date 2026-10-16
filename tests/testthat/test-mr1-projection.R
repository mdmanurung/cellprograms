# MR-1: project_programs_to_cells (replay + projection + canonicalization).

make_mr1_fit <- function(sim, meta = NULL, ...) {
  x <- as_cell_program_data(sim$pseudobulk, sample_metadata = meta)
  fit <- fit_celltype_programs(x, feature_scale = "simulated", ...)
  canonicalize_programs(fit)
}

test_that("projecting the training pseudobulk reproduces canonical donor scores", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  fit <- make_mr1_fit(sim, max_factors = 3, features = "variable", n_variable_genes = 150)
  ct <- fit$cell_types[1]
  cells <- setNames(rep(ct, 24), rownames(fit$fits[[ct]]$Y_used))
  pr <- project_programs_to_cells(fit, sim$pseudobulk[[ct]], cell_types = cells)
  expect_equal(nrow(pr$activity), 24)
  K_all <- sum(vapply(fit$cell_types, function(z) ncol(fit$scores[[z]]), numeric(1)))
  expect_equal(ncol(pr$activity), K_all)   # all programs; off-ct entries NA
  cols <- paste0(ct, "_", seq_len(ncol(fit$scores[[ct]])))
  expect_true(all(cols %in% colnames(pr$activity)))
  for (k in seq_along(cols)) {
    expect_gt(abs(cor(pr$activity[, cols[k]], fit$scores[[ct]][, k])), 0.99)
  }
  expect_equal(pr$coverage$coverage[pr$coverage$cell_type == ct], 1)
  expect_true(all(is.finite(pr$activity[, cols])))   # off-ct columns are NA by design
})

test_that("wrapper equals manual replay + direct oracle + manual canonicalization", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  fit <- make_mr1_fit(sim, max_factors = 3, features = "variable", n_variable_genes = 150)
  ct <- fit$cell_types[1]
  cells <- setNames(rep(ct, 24), rownames(fit$fits[[ct]]$Y_used))
  pr <- project_programs_to_cells(fit, sim$pseudobulk[[ct]], cell_types = cells)
  replay <- fit$fits[[ct]]$preproc_replay
  Q <- as.matrix(sim$pseudobulk[[ct]])[, fit$fits[[ct]]$genes, drop = FALSE]
  Q <- sweep(Q, 2, replay$center, "-")
  Q <- sweep(Q, 2, replay$scale, "/")
  pr2 <- flashier.utils::project_factor_activity(
    fit$fits[[ct]]$flash, query = Q, sample_side = "rows",
    feature_scale = "simulated",
    preprocessing = list(scope = "cellprograms_projection_replay"),
    factors = NULL, ebnm_fn = NULL, maxiter = 200L,
    feature_match = "intersection", min_coverage = 0.9, var_type = 0L,
    dense_budget_bytes = 256 * 1024^2)
  A2 <- sweep(sweep(as.matrix(pr2$representation$activity), 2,
                    fit$canonicalization[[ct]]$scale_multiplier, "/"),
              2, fit$canonicalization[[ct]]$sign_multiplier, "*")
  expect_equal(unname(pr$per_type[[ct]]$activity), unname(A2))
})

test_that("noisy single-cell copies of a donor average back to its score", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  fit <- make_mr1_fit(sim, max_factors = 3, features = "variable", n_variable_genes = 150)
  ct <- fit$cell_types[1]
  Y <- sim$pseudobulk[[ct]]
  nper <- 5L
  expr <- do.call(rbind, lapply(seq_len(nrow(Y)), function(i) {
    noise <- matrix(rnorm(nper * ncol(Y), 0, 0.25 * stats::sd(Y[i, ])), nper, ncol(Y))
    sweep(noise, 2, Y[i, ], "+")
  }))
  rownames(expr) <- paste0("cell_", seq_len(nrow(expr)))
  colnames(expr) <- colnames(Y)
  cells <- setNames(rep(ct, nrow(expr)), rownames(expr))
  pr <- project_programs_to_cells(fit, expr, cell_types = cells)
  donor_mean <- vapply(seq_len(nrow(Y)), function(i) {
    mean(pr$activity[((i - 1L) * nper + 1L):(i * nper), 1])
  }, numeric(1))
  expect_gt(abs(cor(donor_mean, fit$scores[[ct]][, 1])), 0.9)
})

test_that("covariate residualization replays with donors + metadata and refuses without them", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  meta <- data.frame(observation_id = paste0("obs", seq_len(24)),
                     sex = factor(rep(c("F", "M"), each = 12)))
  fit <- make_mr1_fit(sim, meta, max_factors = 3, features = "variable",
                      n_variable_genes = 150, covariates = ~sex)
  ct <- fit$cell_types[1]
  cells <- setNames(rep(ct, 24), rownames(fit$fits[[ct]]$Y_used))
  donors <- setNames(rownames(fit$fits[[ct]]$Y_used), rownames(fit$fits[[ct]]$Y_used))
  pr <- project_programs_to_cells(fit, sim$pseudobulk[[ct]], cell_types = cells,
                                  donors = donors, sample_metadata = meta)
  expect_gt(abs(cor(pr$activity[, 1], fit$scores[[ct]][, 1])), 0.99)
  expect_error(project_programs_to_cells(fit, sim$pseudobulk[[ct]],
    cell_types = cells, donors = donors), "sample_metadata")
  expect_error(project_programs_to_cells(fit, sim$pseudobulk[[ct]],
    cell_types = cells), "donors")
})

test_that("replayed covariate design mismatch fails closed", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  meta <- data.frame(observation_id = paste0("obs", seq_len(24)),
                     sex = factor(rep(c("F", "M"), each = 12)),
                     site = factor(rep("s1", 24), levels = c("s1", "s2")))
  fit <- make_mr1_fit(sim, meta, max_factors = 3, features = "variable",
                      n_variable_genes = 150, covariates = ~sex + site)
  ct <- fit$cell_types[1]
  cells <- setNames(rep(ct, 24), rownames(fit$fits[[ct]]$Y_used))
  donors <- setNames(rownames(fit$fits[[ct]]$Y_used), rownames(fit$fits[[ct]]$Y_used))
  # fit time: constant `site` dropped for this cell type; replay rebuilds the
  # full formula and gets an extra (zero) column -> fail closed
  expect_error(project_programs_to_cells(fit, sim$pseudobulk[[ct]], cell_types = cells,
    donors = donors, sample_metadata = meta), "does not match the design")
})

test_that("projection refuses undeclared units, incompatible units, and fixed-mode fits", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  meta <- data.frame(observation_id = paste0("obs", seq_len(24)),
                     sex = factor(rep(c("F", "M"), each = 12)))
  x_u <- as_cell_program_data(sim$pseudobulk)
  fit_u <- canonicalize_programs(fit_celltype_programs(
    x_u, max_factors = 3, features = "variable", n_variable_genes = 150))
  ct <- fit_u$cell_types[1]
  cells <- setNames(rep(ct, 24), rownames(fit_u$fits[[ct]]$Y_used))
  expect_error(project_programs_to_cells(fit_u, sim$pseudobulk[[ct]],
    cell_types = cells), "undeclared")
  expect_error(project_programs_to_cells(fit_u, sim$pseudobulk[[ct]],
    cell_types = cells, feature_scale = "simulated"), "incompatible")
  fit_fx <- make_mr1_fit(sim, meta, max_factors = 3, features = "variable",
                         n_variable_genes = 150, covariates = ~sex,
                         covariate_mode = "fixed")
  expect_error(project_programs_to_cells(fit_fx, sim$pseudobulk[[ct]],
    cell_types = cells), "cannot be replayed")
  # coverage failure: too few fit genes present in the query
  fit_d <- make_mr1_fit(sim, max_factors = 3, features = "variable", n_variable_genes = 150)
  expr_sub <- sim$pseudobulk[[ct]][, seq_len(20), drop = FALSE]
  expect_error(project_programs_to_cells(fit_d, expr_sub,
    cell_types = setNames(rep(ct, 24), rownames(expr_sub))), "Insufficient|coverage")
})

test_that("K=0 cell types are skipped with NA activities", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  pb <- sim$pseudobulk
  zero <- matrix(5, nrow = 24, ncol = 300,
                 dimnames = list(paste0("z", seq_len(24)), colnames(pb$B)))
  pb$ZERO <- zero
  x <- as_cell_program_data(pb)
  fit <- canonicalize_programs(fit_celltype_programs(
    x, max_factors = 3, features = "variable", n_variable_genes = 150,
    feature_scale = "simulated"))
  expect_true(isTRUE(fit$fits$ZERO$skipped))
  cells <- setNames(c(rep("B", 24), rep("ZERO", 24)),
                    c(rownames(pb$B), rownames(zero)))
  pr <- project_programs_to_cells(fit, rbind(pb$B, zero), cell_types = cells)
  expect_equal(sum(grepl("^ZERO_", colnames(pr$activity))), 0)
  expect_true(all(is.na(pr$activity[cells == "ZERO", , drop = FALSE])))
  expect_equal(pr$coverage$programs[pr$coverage$cell_type == "ZERO"], 0)
})
