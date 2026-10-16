# MR-0: frozen preprocessing specifications (freeze_preprocessing).

test_that("freeze_preprocessing -> spec-governed refit is bit-identical to explicit-argument refit", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150)
  expect_s3_class(spec, "cell_program_preproc_spec")
  expect_true(is.character(spec$digest) && nzchar(spec$digest))
  fit_spec <- fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec)
  fit_exp <- fit_celltype_programs(x, max_factors = 3, features = "variable",
                                   n_variable_genes = 150)
  for (ct in x$cell_types) {
    expect_identical(fit_spec$fits[[ct]]$flash$L_pm, fit_exp$fits[[ct]]$flash$L_pm)
    expect_identical(fit_spec$fits[[ct]]$flash$F_pm, fit_exp$fits[[ct]]$flash$F_pm)
    expect_identical(fit_spec$preprocessing$feature_sets,
                     fit_exp$preprocessing$feature_sets)
  }
  expect_equal(fit_spec$preprocessing$frozen_spec, spec$digest)
  expect_null(fit_exp$preprocessing$frozen_spec)
})

test_that("spec digest is stable across freezes and survives an RDS round-trip", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  s1 <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150)
  s2 <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150)
  expect_identical(s1$digest, s2$digest)
  tmp <- file.path(tempdir(), "mr0-spec.rds")
  saveRDS(s1, tmp)
  s3 <- readRDS(tmp)
  expect_s3_class(s3, "cell_program_preproc_spec")
  expect_identical(s3$digest, s1$digest)
  fit <- fit_celltype_programs(x, max_factors = 3, preprocessing_spec = s3)
  expect_equal(fit$preprocessing$frozen_spec, s1$digest)
})

test_that("spec-governed refit refuses conflicting explicit arguments", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150)
  expect_error(fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec,
                                     features = "variable"), "frozen")
  expect_error(fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec,
                                     n_variable_genes = 100), "frozen")
  expect_error(fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec,
                                     center = FALSE), "frozen")
  expect_error(fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec,
                                     scale = TRUE), "frozen")
  # defaults are not conflicts
  expect_error(fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec), NA)
  # non-spec objects are refused
  expect_error(fit_celltype_programs(x, max_factors = 3, preprocessing_spec = list()),
               "freeze_preprocessing")
})

test_that("frozen covariate formula replays: spec refit == explicit formula refit", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  meta <- data.frame(observation_id = paste0("obs", seq_len(24)),
                     sex = factor(rep(c("F", "M"), each = 12)))
  x <- as_cell_program_data(sim$pseudobulk, sample_metadata = meta)
  spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150,
                               covariates = ~sex)
  expect_equal(spec$covariates, "~sex")
  fit_spec <- fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec)
  fit_exp <- fit_celltype_programs(x, max_factors = 3, features = "variable",
                                   n_variable_genes = 150, covariates = ~sex)
  for (ct in x$cell_types) {
    expect_identical(fit_spec$fits[[ct]]$flash$L_pm, fit_exp$fits[[ct]]$flash$L_pm)
    rr <- fit_spec$fits[[ct]]$preproc_replay
    expect_false(is.null(rr$residualize_coefficients))
    expect_equal(colnames(rr$residualize_coefficients), fit_spec$fits[[ct]]$genes)
    expect_equal(rr$design_columns, colnames(.cov_design(
      rownames(x$pseudobulk[[ct]]), meta, ~sex, ct)))
  }
})

test_that("spec refit fails closed on missing genes or absent cell types", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150)
  sim2 <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 5)
  drop_gene <- spec$feature_sets$B[1]
  x_missing <- as_cell_program_data(lapply(sim2$pseudobulk, function(m) {
    m[, colnames(m) != drop_gene, drop = FALSE]
  }))
  expect_error(fit_celltype_programs(x_missing, max_factors = 3,
                                     preprocessing_spec = spec), "absent from this dataset")
  sim3 <- simulate_pseudobulk(scenario = "private", n_samples = 24, n_genes = 300,
                              program_size = 30, cell_types = c("B", "MAIT"), seed = 2)
  x_other <- as_cell_program_data(sim3$pseudobulk)
  expect_error(fit_celltype_programs(x_other, max_factors = 3,
                                     preprocessing_spec = spec), "lacks cell type")
})

test_that("exclude_artifacts is resolved into the spec and .spec_covariates rehydrates formulas", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150,
                               exclude_artifacts = TRUE)
  expect_true(spec$exclude_artifacts)
  expect_equal(spec$gene_universe, sort(unique(unlist(lapply(x$pseudobulk, colnames)))))
  f <- .spec_covariates("~ age + sex")
  expect_true(inherits(f, "formula"))
  expect_equal(deparse(f), deparse(~ age + sex))
  expect_identical(.spec_covariates(c("age", "sex")), c("age", "sex"))
  expect_null(.spec_covariates(NULL))
})
