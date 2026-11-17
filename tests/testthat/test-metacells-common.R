# Metacell tests: shared helpers of the SuperCell 2.0 + mcRigor path (mcRigor
# output normalization, gamma validation, group ids). The SuperCell builder
# itself is tested in test-metacells-supercell.R.

# ---- fixtures ---------------------------------------------------------------
fake_mcrigor_res <- function(n_mc = 5) {
  # pinned against mcRigor_DETECT's return (ledger 2026-10-18, bridge fix 6):
  # TabMC with T_org + mcRigor verdict column; thre data.frame with metacell
  # ids as rownames.
  tab <- data.frame(
    gamma = rep(25, n_mc),
    ZeroRate = runif(n_mc, 0, 0.3),
    size = sample(40:120, n_mc),
    T_org = runif(n_mc, 0, 0.5),
    T_colperm = runif(n_mc, 0, 0.5),
    TT_div = runif(n_mc, 0, 0.5),
    mcRigor = sample(c("trustworthy", "dubious"), n_mc, replace = TRUE),
    stringsAsFactors = FALSE
  )
  rownames(tab) <- paste0("g", seq_len(n_mc))
  thre <- data.frame(size = tab$size, thre = runif(n_mc, 0.3, 0.6),
                     stringsAsFactors = FALSE)
  rownames(thre) <- paste0("g", seq_len(n_mc))
  list(TabMC = tab, thre = thre)
}

# ---- .mcrigor_normalize ------------------------------------------------------
test_that("mcRigor_DETECT output normalizes to per-metacell QC table", {
  res <- fake_mcrigor_res(6)
  out <- cellprograms:::.mcrigor_normalize(res, group_id = "donor1_CD8")
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 6L)
  expect_setequal(colnames(out),
                  c("metacell", "Group", "mcDiv", "mcDiv_threshold", "dubious"))
  expect_true(all(out$Group == "donor1_CD8"))
  expect_identical(out$metacell, paste0("g", 1:6))
  expect_true(all(is.finite(out$mcDiv)))
  expect_true(all(is.finite(out$mcDiv_threshold)))
  expect_identical(out$dubious, res$TabMC$mcRigor != "trustworthy")
})

test_that("dubious falls back to mcDiv > threshold when no verdict column", {
  res <- fake_mcrigor_res(5)
  res$TabMC$mcRigor <- NULL
  res$TabMC$T_org <- c(0.9, 0.1, 0.5, 0.4, 0.2)
  out <- cellprograms:::.mcrigor_normalize(res, group_id = "g")
  expect_identical(out$dubious, out$mcDiv > out$mcDiv_threshold)
})

test_that("unrecognized mcRigor return structure raises loudly", {
  expect_error(cellprograms:::.mcrigor_normalize(list(a = 1), group_id = "g"),
               "unrecognized")
})

test_that(".mc_validate_gammas dedups, sorts, and rejects bad input", {
  expect_identical(.mc_validate_gammas(c(10, 5, 10)), c(5L, 10L))
  expect_identical(.mc_validate_gammas(7), 7L)
  expect_error(.mc_validate_gammas("a"), "positive integers")
  expect_error(.mc_validate_gammas(numeric(0)), "positive integers")
  expect_error(.mc_validate_gammas(c(5, NA)), "positive integers")
  expect_error(.mc_validate_gammas(c(5, -1)), "positive integers")
  expect_error(.mc_validate_gammas(c(5, 2.5)), "positive integers")
})

test_that(".mc_get_group_ids fails closed on NA or missing group.by columns", {
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  counts <- Matrix::Matrix(rpois(40, 2), nrow = 10, ncol = 4, sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:10)
  colnames(counts) <- paste0("cell", 1:4)
  meta <- data.frame(donor = c("d1", "d1", NA, "d2"), row.names = colnames(counts))
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  expect_error(.mc_get_group_ids(obj, "donor"), "NA in group.by")
  expect_error(.mc_get_group_ids(obj, "no_such_col"), "not found")
  # clean metadata still works
  meta$donor[is.na(meta$donor)] <- "d2"
  obj2 <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  ids <- .mc_get_group_ids(obj2, "donor")
  expect_identical(unname(ids), c("d1", "d1", "d2", "d2"))
})
