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

# ---- OptimizeGammaSEACells membership coding (Seurat v5 mangling guard) ----
test_that("membership coding drops single-metacell gammas and integer-codes the rest", {
  cells <- paste0("cell", 1:50)
  mklab <- function(sizes) rep(paste0("mc_", seq_along(sizes)), sizes)
  mem_cols <- list(
    mklab(c(13, 12, 13, 12)),   # gamma 25: 4 metacells
    mklab(c(25, 25)),           # gamma 50: 2 metacells
    mklab(c(50))                # gamma 100: collapsed to 1 -> not scoreable
  )
  coded <- cellprograms:::.mc_optimize_membership_codes(mem_cols, c(25, 50, 100), cells)
  expect_equal(coded$keep_gamma, c(25, 50))
  expect_equal(coded$drop_gamma, 100)
  expect_equal(colnames(coded$code_df), c("25", "50"))
  expect_equal(rownames(coded$code_df), cells)
  # codes are per-column integers 1..K (buildmc's g-strip contract)
  expect_equal(sort(unique(coded$code_df[["25"]])), 1:4)
  expect_equal(sort(unique(coded$code_df[["50"]])), 1:2)
  # labels are recoverable: code 1 maps back to mc_1 in every column
  expect_equal(unique(mem_cols[[1]][coded$code_df[["25"]] == 1]), "mc_1")
  # all gammas collapsed -> nothing to optimize
  collapsed <- list(mklab(c(50)), mklab(c(50)))
  coded0 <- cellprograms:::.mc_optimize_membership_codes(collapsed, c(50, 100), cells)
  expect_equal(length(coded0$keep_gamma), 0L)
  expect_null(coded0$code_df)
})

test_that("GetMetacellSEACells returns stored results or stops informatively", {
  # 2 genes x 1 cell: a 1x1 sparse matrix is diagonal (ddiMatrix) and Seurat's
  # coercion warning would pollute the suite warning tally
  obj <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(c(1, 2), 2, 1, sparse = TRUE,
                            dimnames = list(c("g1", "g2"), "c1")))
  expect_error(GetMetacellSEACells(obj), "run MetacellsByGroupsSEACells first")
  obj@misc$SEACells <- list(metacell_object = "x", qc = data.frame())
  out <- GetMetacellSEACells(obj)
  expect_identical(out$metacell_object, "x")
  # custom storage name
  obj@misc$other <- list(y = 2)
  expect_identical(GetMetacellSEACells(obj, wgcna_name = "other")$y, 2)
})

test_that("mcRigor_OPTIMIZE accepts integer-coded membership (Seurat v5 regression)", {
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  skip_if_not(requireNamespace("mcRigor", quietly = TRUE))
  set.seed(11)
  counts <- Matrix::Matrix(rpois(360, 2), nrow = 6, ncol = 60, sparse = TRUE,
                           dimnames = list(paste0("gene", 1:6), paste0("cell", 1:60)))
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  # Integer codes, the post-fix contract: buildmc aligns sizes through
  # AggregateExpression colnames ("g1,g2,g3" -> strip -> "1","2","3").
  # Pre-fix, "mc_<k>" labels mangled to "mc-1,..." died in buildmc with
  # "No cell overlap between new meta data and Seurat object".
  lab <- rep(1:3, each = 20)
  code_df <- data.frame(`5` = lab, `10` = lab, check.names = FALSE,
                        stringsAsFactors = FALSE)
  # buildmc derives cell IDs from membership rownames (the wrapper sets
  # rownames(df) <- cells_g); without them the cell subset is empty
  rownames(code_df) <- paste0("cell", 1:60)
  args <- list(Nrep = 5, gene_filter = 0.1, feature_use = 2000,
               test_cutoff = 0.01, cor_method = "pearson",
               aggregate_method = "mean", draw = FALSE, step_save = FALSE,
               output_file = tempfile(fileext = ".rds"), fields = character(0))
  # mcRigor emits numerical chatter on random data (documented in the ledger);
  # the assertion is on the result structure.
  res <- suppressWarnings(do.call(mcRigor::mcRigor_OPTIMIZE,
    utils::modifyList(list(obj_singlecell = obj, cell_membership = code_df), args)))
  best <- cellprograms:::.mcrigor_optimize_pick(res, c(5, 10))
  expect_length(best$scores, 2L)
  expect_true(is.na(best$gamma) || best$gamma %in% c(5, 10))
})

# Planted-archetype fixture: 5 subpopulations per 50-cell stratum, each with a
# private 6-gene marker block (lambda 8 vs 0.5 background), so SEACells has
# real structure to find and the mcRigor QC path runs on genuine multi-metacell
# output. The default cells_per_metacell = 75 gives 1 metacell per group and
# skips QC entirely (the single-metacell edge case); random data collapses to
# 1 metacell even at cells_per_metacell = 10 (see the collapse test below).
.planted_metacell_fixture <- function() {
  set.seed(5)
  lam <- matrix(0.5, nrow = 60, ncol = 50)
  for (s in 1:5) {
    cells <- ((s - 1L) * 10L + 1L):(s * 10L)
    genes <- ((s - 1L) * 6L + 1L):(s * 6L)
    lam[genes, cells] <- 8
  }
  counts <- Matrix::Matrix(
    rpois(60 * 200, as.vector(do.call(cbind, rep(list(lam), 4)))),
    nrow = 60, ncol = 200, sparse = TRUE
  )
  rownames(counts) <- paste0("gene", 1:60)
  colnames(counts) <- paste0("cell", 1:200)
  meta <- data.frame(
    donor = rep(c("d1", "d2"), each = 100),
    group = rep(c("A", "B"), times = 100)
  )
  rownames(meta) <- colnames(counts)
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 5, verbose = FALSE)
  obj
}

test_that("integration: mcRigor QC and gamma optimization on planted multi-metacell groups", {
  skip_if_not(identical(Sys.getenv("CPMC_INTEGRATION"), "1"),
              "set CPMC_INTEGRATION=1 to run the basilisk/SEACells integration test")
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  obj <- .planted_metacell_fixture()

  # QC path: mcRigor_DETECT runs on real SEACells output with >= 2 metacells
  out <- MetacellsByGroupsSEACells(obj, group.by = c("donor", "group"),
                                   n_top_genes = 30, n_pc = 5,
                                   min_cells = 50, cells_per_metacell = 10,
                                   seed = 5, verbose = FALSE)
  mc <- out@misc$SEACells
  expect_equal(nrow(mc$skipped_groups), 0L)
  expect_setequal(unique(mc$qc$Group), c("d1_A", "d1_B", "d2_A", "d2_B"))
  expect_true(all(table(mc$qc$Group) >= 2L))  # fixture point: multi-metacell groups
  expect_true(any(is.finite(mc$qc$mcDiv)))  # mcRigor ran (fallback rows are NA)

  # Optimization path: 2 donor groups x 2 gammas; pre-fix, every group died in
  # buildmc's size assignment under Seurat v5.
  opt <- OptimizeGammaSEACells(out, group.by = "donor", gammas = c(5, 10),
                               n_top_genes = 30, n_pc = 5, min_cells = 50,
                               mcrigor_args = list(Nrep = 5), seed = 5,
                               verbose = FALSE)
  expect_true(length(opt$optimal_gamma) > 0)
  expect_true(all(opt$optimal_gamma %in% c(5, 10)))
  expect_true(nrow(opt$scores) >= 2L)
})

test_that("integration: archetype collapse drops misaligned soft assignments gracefully", {
  skip_if_not(identical(Sys.getenv("CPMC_INTEGRATION"), "1"),
              "set CPMC_INTEGRATION=1 to run the basilisk/SEACells integration test")
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  set.seed(5)
  counts <- Matrix::Matrix(rpois(600, 2), nrow = 60, ncol = 100, sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:60)
  colnames(counts) <- paste0("cell", 1:100)
  meta <- data.frame(donor = rep("d1", 100))
  rownames(meta) <- colnames(counts)
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 30, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 5, verbose = FALSE)
  # Random data: SEACells fits 10 archetypes but the hard assignments collapse
  # (diagnostic: all cells -> 1 metacell). Pre-fix this crashed naming the
  # archetype-columned soft matrix with metacell names.
  out <- suppressWarnings(MetacellsByGroupsSEACells(
    obj, group.by = "donor", n_top_genes = 30, n_pc = 5, min_cells = 50,
    cells_per_metacell = 10, seed = 5, verbose = FALSE))
  mc <- out@misc$SEACells
  expect_equal(nrow(mc$skipped_groups), 0L)
  expect_true(nrow(mc$qc) < 10L)  # collapse happened (10 archetypes requested)
  expect_true(is.null(mc$soft_assignments))  # misaligned soft matrix dropped
})

# ---- edge cases (2026-10-20 segment) ----------------------------------------

test_that(".mc_validate_gammas dedups, sorts, and rejects bad input", {
  expect_identical(.mc_validate_gammas(c(10, 5, 10)), c(5L, 10L))
  expect_identical(.mc_validate_gammas(7), 7L)
  expect_error(.mc_validate_gammas("a"), "positive integers")
  expect_error(.mc_validate_gammas(numeric(0)), "positive integers")
  expect_error(.mc_validate_gammas(c(5, NA)), "positive integers")
  expect_error(.mc_validate_gammas(c(5, -1)), "positive integers")
  expect_error(.mc_validate_gammas(c(5, 2.5)), "positive integers")
})

test_that("membership coding excludes NA-labeled cells and never counts NA as a level", {
  cells <- paste0("c", 1:6)
  mem <- list(
    "5"  = c("a", "a", "b", "b", NA, NA),   # 2 real labels + NA cells
    "10" = c("a", "a", "a", "a", NA, "a")   # 1 real label once NA excluded
  )
  expect_warning(
    coded <- .mc_optimize_membership_codes(mem, c(5, 10), cells),
    "NA membership labels excluded"
  )
  # gamma 10 is degenerate once NA cells are gone (unique() would have
  # counted NA as a second level pre-fix)
  expect_identical(coded$keep_gamma, 5)
  expect_identical(coded$drop_gamma, 10)
  expect_identical(rownames(coded$code_df), cells[1:4])
  expect_false(any(is.na(coded$code_df[[1]])))
})

test_that("membership coding keeps a single valid gamma and preserves input order", {
  cells <- paste0("c", 1:4)
  mem <- list("10" = c("a", "a", "b", "b"), "5" = c("x", "x", "y", "y"))
  coded <- .mc_optimize_membership_codes(mem, c(10, 5), cells)
  expect_identical(coded$keep_gamma, c(10, 5))   # input order, not sorted
  expect_identical(colnames(coded$code_df), c("10", "5"))
  single <- .mc_optimize_membership_codes(mem[1], 10, cells)
  expect_identical(single$keep_gamma, 10)
  expect_identical(sort(unique(single$code_df[[1]])), c(1L, 2L))
})

test_that("membership coding with all-degenerate gammas returns NULL code_df", {
  cells <- paste0("c", 1:4)
  mem <- list("5" = rep("a", 4), "10" = rep("b", 4))
  coded <- .mc_optimize_membership_codes(mem, c(5, 10), cells)
  expect_identical(coded$keep_gamma, numeric(0))
  expect_identical(coded$drop_gamma, c(5, 10))
  expect_null(coded$code_df)
})

test_that(".mcrigor_optimize_pick is honest about all-NA scores and single gamma", {
  res <- data.frame(gamma = c(5, 10), Score = c(NA_real_, NA_real_))
  pick <- .mcrigor_optimize_pick(res, c(5, 10))
  expect_true(is.na(pick$gamma))
  expect_true(all(is.na(pick$scores)))
  one <- .mcrigor_optimize_pick(data.frame(gamma = 5, Score = 0.4), 5)
  expect_identical(one$gamma, 5)
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

test_that("integration: single-gamma optimization and no-prior-run optimize", {
  skip_if_not(identical(Sys.getenv("CPMC_INTEGRATION"), "1"),
              "set CPMC_INTEGRATION=1 to run the basilisk/SEACells integration test")
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  obj <- .planted_metacell_fixture()
  # single candidate gamma: the pick must be that gamma for every group
  opt1 <- OptimizeGammaSEACells(obj, group.by = "donor", gammas = 10,
                                n_top_genes = 30, n_pc = 5, min_cells = 50,
                                mcrigor_args = list(Nrep = 5), seed = 5,
                                verbose = FALSE)
  expect_true(all(unlist(opt1$optimal_gamma) == 10))
  # OptimizeGammaSEACells refits SEACells per gamma and does not require a
  # prior MetacellsByGroupsSEACells run on the object
  fresh <- .planted_metacell_fixture()
  expect_null(fresh@misc$SEACells)
  opt2 <- OptimizeGammaSEACells(fresh, group.by = "donor", gammas = c(5, 10),
                                n_top_genes = 30, n_pc = 5, min_cells = 50,
                                mcrigor_args = list(Nrep = 5), seed = 5,
                                verbose = FALSE)
  expect_true(length(opt2$optimal_gamma) > 0)
})

test_that("integration: all groups below min_cells aborts fail-closed", {
  skip_if_not(identical(Sys.getenv("CPMC_INTEGRATION"), "1"),
              "set CPMC_INTEGRATION=1 to run the basilisk/SEACells integration test")
  skip_if_not(requireNamespace("Seurat", quietly = TRUE))
  obj <- .planted_metacell_fixture()
  expect_error(
    MetacellsByGroupsSEACells(obj, group.by = c("donor", "group"),
                              n_top_genes = 30, n_pc = 5,
                              min_cells = 10000, cells_per_metacell = 10,
                              seed = 5, verbose = FALSE),
    "no groups passed min_cells"
  )
})
