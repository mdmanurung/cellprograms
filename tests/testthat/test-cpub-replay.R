test_that("CPUB-002 trained designs are subset/repeat/reorder invariant", {
  meta <- data.frame(observation_id = paste0("d", 1:40),
    age = seq(20, 80, length.out = 40), sex = factor(rep(c("F", "M"), 20)))
  forms <- list(~age + sex, ~poly(age, 2) + sex, ~scale(age) + sex,
                ~splines::ns(age, df = 3) + sex, ~splines::bs(age, df = 4),
                ~log(age), ~poly(age, 2, raw = TRUE))
  for (f in forms) {
    X <- .cov_design(meta$observation_id, meta, f, "A")
    state <- attr(X, "prediction_state")
    for (idx in list(c(7, 2, 7), 12L, 40:1, 1:20 * 2L)) {
      p <- .replay_design(meta$observation_id[idx], meta, f, colnames(X), "A", state)
      expect_equal(unname(p), unname(X[idx, , drop = FALSE]),
                   ignore_attr = TRUE, tolerance = 1e-12)
    }
    copy <- unserialize(serialize(state, NULL))
    p <- .predict_covariate_design(copy, meta$observation_id[1:3], meta, "A")
    expect_equal(unname(p), unname(X[1:3, , drop = FALSE]), ignore_attr = TRUE)
  }
})

test_that("CPUB-002 constant terms stay dropped and new levels fail closed", {
  meta <- data.frame(observation_id = paste0("d", 1:40), age = 1:40,
    site = factor(rep("a", 40)), sex = rep(c("F", "M"), 20))
  X <- .cov_design(meta$observation_id, meta, ~age + site + age:site + sex, "A")
  state <- attr(X, "prediction_state")
  expect_false("site" %in% state$variables)
  p <- .predict_covariate_design(state, rep("d1", 3), meta, "A")
  expect_equal(unname(p), unname(X[rep(1, 3), , drop = FALSE]), ignore_attr = TRUE)
  meta$sex[1] <- "new"
  expect_error(.predict_covariate_design(state, "d1", meta, "A"), "new level")
  expect_error(.replay_design("d1", meta, ~age, colnames(X), "A"), "refit legacy")
  expect_error(.cov_design(meta$observation_id, meta, ~rank(age), "A"),
               "Unsupported covariate transform")
  expect_error(.cov_design(meta$observation_id, meta, ~scale(log(age)), "A"),
               "Only scale")
})

test_that("CPUB-011 frozen specs verify contents, versions and units", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 100, program_size = 10, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  s <- freeze_preprocessing(x, feature_scale = "log1p_cp10k")
  expect_identical(s$spec_version, "1.1.0")
  expect_identical(s$feature_scale, "log1p_cp10k")
  expect_error(.validate_preproc_spec(s), NA)
  for (field in c("center", "scale")) {
    bad <- s; bad[[field]] <- !bad[[field]]
    expect_error(.validate_preproc_spec(bad), "digest mismatch")
  }
  bad <- s; bad$feature_sets[[1]] <- rev(bad$feature_sets[[1]])
  expect_error(.validate_preproc_spec(bad), "digest mismatch")
  bad <- s; bad$spec_version <- "9"
  expect_error(.validate_preproc_spec(bad), "Unsupported")
  # Old intact specs are accepted with unknown units; never silently upgrade.
  old <- s; old$digest <- NULL; old$spec_version <- "1.0.0"
  old$feature_scale <- NULL; old$digest <- digest::digest(old, algo = "sha256")
  expect_error(.validate_preproc_spec(old), NA)
})

test_that("CPUB-002 trained transformed-covariate subtraction reproduces Y_used", {
  set.seed(302)
  ids <- paste0("d", 1:40)
  meta <- data.frame(observation_id = ids, age = seq(20, 80, length.out = 40),
                     sex = factor(rep(c("F", "M"), 20)))
  y <- matrix(rnorm(40 * 30), 40, 30,
              dimnames = list(ids, paste0("g", 1:30)))
  x <- as_cell_program_data(list(A = y), sample_metadata = meta)
  f <- fit_celltype_programs(x, covariates = ~poly(age, 2) + scale(age) + sex,
                             max_factors = 1, feature_scale = "diagnostic")
  rr <- f$fits$A$preproc_replay
  # The alias between scale(age) and the polynomial linear term is deliberate.
  beta <- rr$residualize_coefficients; beta[is.na(beta)] <- 0
  for (idx in list(c(7, 2, 7), 40:1, 3L)) {
    Xq <- .predict_covariate_design(rr$design_state, ids[idx], meta, "A")
    q <- y[idx, f$fits$A$genes, drop = FALSE] -
      Xq %*% beta[rr$design_columns, f$fits$A$genes, drop = FALSE]
    q <- sweep(sweep(q, 2, rr$center, "-"), 2, rr$scale, "/")
    expect_equal(unname(q), unname(f$fits$A$Y_used[idx, , drop = FALSE]),
                 tolerance = 1e-10)
  }
})
