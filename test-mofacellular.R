test_that("as_mofacellular_input produces create_init_exp-ready counts/coldata", {
  set.seed(5)
  n <- 8; g <- 20
  Y1 <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  Y2 <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  cpd <- as_cell_program_data(list(A = Y1, B = Y2),
                              sample_metadata = data.frame(observation_id = paste0("s", 1:n)))
  out <- as_mofacellular_input(cpd)

  # shape: genes rows, profiles (cell types x samples) columns
  expect_equal(dim(out$counts), c(g, 2 * n))
  expect_equal(ncol(out$coldata), 3)
  expect_equal(sort(names(out$coldata)), c("cell_counts", "cell_type", "donor_id"))
  expect_true(all(out$coldata$cell_counts == 8L) || all(is.na(out$coldata$cell_counts)))

  # column names and rownames alignment
  expect_setequal(colnames(out$counts), c(paste0("A_s", 1:n), paste0("B_s", 1:n)))
  expect_equal(rownames(out$coldata), colnames(out$counts))

  # values round-trip: counts block for cell type A equals t(Y1)
  expect_equal(unname(out$counts[, paste0("A_s", 1:n)]), unname(t(Y1)))
  expect_equal(unname(out$counts[, paste0("B_s", 1:n)]), unname(t(Y2)))

  # works on a fit object too (pseudobulk carried through)
  skip_if_not_installed("flashier")
  fit <- fit_celltype_programs(cpd, seed = 1)
  out2 <- as_mofacellular_input(fit)
  expect_equal(out2$counts, out$counts)

  # explicit cell_counts matching
  cc <- rbind(data.frame(observation_id = paste0("s", 1:n), cell_type = "A", cell_counts = 10),
              data.frame(observation_id = paste0("s", 1:n), cell_type = "B", cell_counts = 20))
  out3 <- as_mofacellular_input(cpd, cell_counts = cc)
  expect_equal(sort(out3$coldata$cell_counts), sort(rep(c(10, 20), each = n)))
  expect_error(as_mofacellular_input(cpd, cell_counts = cc[1:5, ]), "missing")
})

test_that("as_mofacellular_input handles missing cell-type observations", {
  set.seed(6)
  g <- 15
  Y1 <- matrix(rnorm(6 * g), 6, g, dimnames = list(paste0("s", 1:6), paste0("g", 1:g)))
  Y2 <- matrix(rnorm(4 * g), 4, g, dimnames = list(paste0("s", 3:6), paste0("g", 1:g)))
  cpd <- as_cell_program_data(list(A = Y1, B = Y2))
  out <- as_mofacellular_input(cpd)
  expect_equal(ncol(out$counts), 6 + 4)
  expect_equal(sum(out$coldata$cell_type == "B"), 4)
})
