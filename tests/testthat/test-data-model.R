test_that("as_cell_program_data accepts a valid multi-cell-type input with missing observations", {
  set.seed(1)
  genes <- paste0("G", 1:50)
  ids_all <- paste0("obs", 1:10)
  mk <- function(n) matrix(rnorm(n * 50), nrow = n, ncol = 50,
                           dimnames = list(ids_all[seq_len(n)], genes))
  pb <- list(B = mk(10), CD4 = mk(10), NK = mk(7))  # NK missing 3 observations
  meta <- data.frame(observation_id = ids_all, participant_id = ids_all)

  x <- as_cell_program_data(pb, sample_metadata = meta)
  expect_s3_class(x, "cell_program_data")
  expect_setequal(x$cell_types, c("B", "CD4", "NK"))
  expect_length(x$shared_genes, 50)
  expect_message(print(x), "NK")
})

test_that("data contract failures fail loudly", {
  genes <- paste0("G", 1:10)
  ids <- paste0("obs", 1:5)
  m <- matrix(rnorm(50), 5, 10, dimnames = list(ids, genes))
  meta <- data.frame(observation_id = ids)

  # duplicate observation IDs within a cell type
  bad <- m; rownames(bad)[2] <- "obs1"
  expect_error(as_cell_program_data(list(B = bad), meta), "Duplicated observation")

  # duplicate genes
  badg <- m; colnames(badg)[2] <- "G1"
  expect_error(as_cell_program_data(list(B = badg), meta), "Duplicated gene")

  # duplicate cell-type names
  expect_error(as_cell_program_data(list(B = m, B = m), meta), "Duplicated cell-type")

  # NA values
  badna <- m; badna[1, 1] <- NA
  expect_error(as_cell_program_data(list(B = badna), meta), "NA")

  # metadata missing observation_id
  expect_error(as_cell_program_data(list(B = m), data.frame(x = 1:5)), "observation_id")

  # observation in matrix but not metadata
  expect_error(as_cell_program_data(list(B = m), data.frame(observation_id = ids[1:4])),
               "missing from")
})

test_that("validate warns on small N and zero-variance genes", {
  genes <- paste0("G", 1:10)
  ids <- paste0("obs", 1:5)
  m <- matrix(rnorm(50), 5, 10, dimnames = list(ids, genes))
  m[, 1] <- 1  # zero variance gene
  x <- as_cell_program_data(list(B = m),
                            data.frame(observation_id = ids))
  expect_warning(validate_cell_program_data(x), "near-zero variance")
  expect_warning(validate_cell_program_data(x), "only 5 observations")
})

test_that("package defaults match the optimized configuration (eval-v3 selection)", {
  # Optimized by the tusoskill benchmark-first loop (run cellprograms-ebmf-config-v1):
  # per-gene variance + point-Laplace gene prior; all other plan Section 9 defaults.
  expect_equal(formals(fit_celltype_programs)$loading_prior, "point_laplace")
  expect_equal(formals(fit_celltype_programs)$var_type, 1)
  expect_equal(formals(fit_celltype_programs)$center, TRUE)
  expect_equal(formals(fit_celltype_programs)$backfit, TRUE)
  expect_equal(formals(fit_celltype_programs)$nullcheck, TRUE)
})
