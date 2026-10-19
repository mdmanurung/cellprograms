# Metacell (SEACells + mcRigor) tests — folded in from the retired
# cellprogramsMC sub-package (ledger 2026-10-18).
#
# R-side logic is tested directly (mcRigor output normalization, gamma-score
# extraction, hard-label reconstruction, export round trips, feasibility
# rule). The full SEACells builder path requires the basilisk environment
# (heavy: env construction + SEACells fit) and is skip-guarded behind
# CPMC_INTEGRATION=1; the smoke-level builder behavior is covered by the W5
# battery and the recorded smoke test (ledger 2026-10-18).

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

# ---- .mcrigor_optimize_pick ---------------------------------------------------
test_that("gamma optimization picks the max-score gamma", {
  gammas <- c(25, 50, 100)
  res <- data.frame(gamma = gammas, Score = c(0.4, 0.9, 0.6))
  best <- cellprograms:::.mcrigor_optimize_pick(res, gammas)
  expect_identical(best$gamma, 50)
  expect_identical(as.numeric(best$scores), c(0.4, 0.9, 0.6))
})

test_that("gamma optimization survives a score-only table", {
  gammas <- c(25, 50, 100)
  res <- data.frame(opt_score = c(0.2, 0.3, 0.8))
  best <- cellprograms:::.mcrigor_optimize_pick(res, gammas)
  expect_identical(best$gamma, 100)
})

# ---- .mc_n_seacells_for_group -------------------------------------------------
test_that("metacell count feasibility rule", {
  expect_identical(cellprograms:::.mc_n_seacells_for_group(1000, 25, NULL), 40L)
  expect_identical(cellprograms:::.mc_n_seacells_for_group(1000, 25, 10), 10L)
  expect_identical(cellprograms:::.mc_n_seacells_for_group(10, 25, NULL), 1L)
})

# ---- .hard_labels_vector ------------------------------------------------------
test_that("bare soft-assignment MTX gets cell names and argmax labels", {
  # bridge fix 8: soft_assignments.mtx is a bare MTX (no dimnames); rows are
  # in exported cell order
  cells <- paste0("c", 1:7)
  soft <- Matrix::Matrix(c(0.9, 0.1, 0.2, 0.8, 0.5, 0.5, 0.1,
                           0.1, 0.9, 0.8, 0.2, 0.5, 0.5, 0.9),
                         nrow = 7, ncol = 2, sparse = TRUE)
  mc_out <- list(soft = soft)
  lab <- cellprograms:::.hard_labels_vector(mc_out, cells, paste0("mc_", 1:2))
  expect_identical(as.character(lab[1]), "mc_1")
  expect_identical(as.character(lab[2]), "mc_2")
  expect_identical(as.character(lab[7]), "mc_2")
  expect_identical(names(lab), cells)  # NA-named membership fixed
})

test_that("missing soft assignments raise an informative error", {
  expect_error(cellprograms:::.hard_labels_vector(list(soft = NULL), letters[1:3],
                                                paste0("mc_", 1:2)),
               "soft assignments required")
})

# ---- export round trip ---------------------------------------------------------
test_that("group export writes the files the python helper consumes", {
  counts <- Matrix::Matrix(rpois(30, 5), nrow = 6, ncol = 5, sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:6)
  colnames(counts) <- paste0("cell", 1:5)
  emb <- matrix(rnorm(10), nrow = 5, dimnames = list(paste0("cell", 1:5), c("PC1", "PC2")))
  paths <- cellprograms:::.export_group_for_python(counts, emb, tempdir())
  expect_true(all(file.exists(unlist(paths))))
  # round-trip fidelity is the contract (the python helper reads via scipy and
  # takes names from cells.csv/genes.csv — MTX carries no dimnames);
  # readMM's sparse-representation class (dgC vs dgT) is Matrix-version detail
  back <- Matrix::readMM(paths$mtx)
  expect_s4_class(back, "sparseMatrix")
  expect_equal(unname(as.matrix(back)), unname(as.matrix(counts)))
})

# ---- integration (heavy; opt-in) ------------------------------------------------
test_that("MetacellsByGroupsSEACells end-to-end on a tiny synthetic object", {
  skip_if_not(identical(Sys.getenv("CPMC_INTEGRATION"), "1"),
              "set CPMC_INTEGRATION=1 to run the basilisk/SEACells integration test")
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  set.seed(5)
  counts <- Matrix::Matrix(rpois(600, 2), nrow = 60, ncol = 200, sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:60)
  colnames(counts) <- paste0("cell", 1:200)
  meta <- data.frame(
    donor = rep(c("d1", "d2"), each = 100),
    group = rep(c("A", "B"), times = 100)
  )
  rownames(meta) <- colnames(counts)
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  # Seurat v5's PrepDR5 errors without variable features (older versions
  # silently used all genes); select them explicitly.
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 5, verbose = FALSE)
  out <- MetacellsByGroupsSEACells(obj, group.by = c("donor", "group"),
                                   n_top_genes = 30, n_pc = 5,
                                   min_cells = 50, seed = 5, verbose = FALSE)
  expect_false(is.null(out@misc$SEACells))
  # skipped_groups is a data.frame: count ROWS (length() would count columns)
  expect_true(nrow(out@misc$SEACells$skipped_groups) < 2)
})
