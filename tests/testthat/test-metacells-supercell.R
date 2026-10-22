# Tests for the SuperCell 2.0 metacell backend (pure R — runs in the BASELINE
# suite; no basilisk/Python needed, unlike the SEACells e2e).
#
# The mcRigor calling convention is pinned against mcRigor HEAD c70cd6f by the
# benchmark (benchmarks/metacell_comparison/): DETECT requires MC-prefixed
# STRING labels named by cell id (integer labels coerce to NA in DETECT's
# named_membership filter -> "Cannot find cells provided"); OPTIMIZE accepts
# both integer codes and MC strings.

.sc_test_fixture <- function(n_cells = 200, n_genes = 600) {
    set.seed(5)
    counts <- matrix(rpois(n_genes * n_cells, 2), nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("gene", seq_len(n_genes))
    colnames(counts) <- paste0("cell", seq_len(n_cells))
    counts[1:20, 1:(n_cells / 2)] <- rpois(20 * (n_cells / 2), 8)
    counts[21:40, ((n_cells / 2) + 1):n_cells] <- rpois(20 * (n_cells / 2), 8)
    counts <- Matrix::Matrix(counts, sparse = TRUE)
    meta <- data.frame(
        donor = rep(c("d1", "d2"), each = n_cells / 2),
        group = rep(c("A", "B"), times = n_cells / 2)
    )
    rownames(meta) <- colnames(counts)
    obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 100, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = 5, verbose = FALSE)
    obj
}

test_that("SCimplify gamma bisection hits the target metacell count", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    obj <- .sc_test_fixture()
    obj_g <- subset(obj, cells = colnames(obj)[1:100])
    hit <- cellprograms:::.sc_gamma_for_target(
        obj_g, target = 10, reduction = "pca", dims = 1:5, assay = "RNA",
        k.knn = 30
    )
    expect_false(is.null(hit$mem))
    expect_lte(abs(hit$n - 10), 0.10 * 10)
    # higher target -> more metacells
    hit2 <- cellprograms:::.sc_gamma_for_target(
        obj_g, target = 25, reduction = "pca", dims = 1:5, assay = "RNA",
        k.knn = 30
    )
    expect_gt(hit2$n, hit$n)
})

test_that("metacell count aggregation: sum and average modes", {
    # byrow = TRUE: Matrix fills column-major by default; the hand-computed
    # expectations below assume g1 = (1, 2, 3), g2 = (4, 5, 6).
    counts <- Matrix::Matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE,
                             sparse = TRUE,
                             dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
    mem <- c(1, 1, 2); names(mem) <- colnames(counts)
    sum_mc <- cellprograms:::.sc_aggregate_counts(counts, mem, c("m1", "m2"), "sum")
    expect_equal(as.numeric(sum_mc["g1", ]), c(3, 3))   # 1+2, 3
    expect_equal(as.numeric(sum_mc["g2", ]), c(9, 6))   # 4+5, 6
    avg_mc <- cellprograms:::.sc_aggregate_counts(counts, mem, c("m1", "m2"), "average")
    expect_equal(as.numeric(avg_mc["g1", ]), c(1.5, 3)) # (1+2)/2, 3/1
    expect_equal(as.numeric(avg_mc["g2", ]), c(4.5, 6))
})

test_that("MetacellsByGroupsSuperCell end-to-end on a tiny synthetic object", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    obj <- .sc_test_fixture()
    out <- MetacellsByGroupsSuperCell(
        obj, group.by = c("donor", "group"), min_cells = 50,
        cells_per_metacell = 10, seed = 5, verbose = FALSE
    )
    res <- out@misc$SuperCell
    expect_false(is.null(res))
    expect_true(is.null(res$soft_assignments))   # hard-membership backend
    # both donor x group combinations processed (100 cells each >= min_cells)
    expect_equal(nrow(res$skipped_groups), 0L)
    # metacell object carries Group + n_cells + the individual group columns
    mo <- res$metacell_object
    expect_true(!is.null(mo[["Group"]]))
    expect_true(!is.null(mo[["n_cells"]]))
    expect_true(all(c("donor", "group") %in% colnames(mo@meta.data)))
    # metacell names are <group>_<i> per group
    expect_true(all(grepl("^donor_group_", colnames(mo))))
    # qc table aligned with the metacell object
    expect_equal(nrow(res$qc), ncol(mo))
    expect_true(all(res$qc$metacell %in% colnames(mo)))
    # fit_stats records gamma and metacell count per group
    expect_true(all(c("gamma", "n_metacells") %in% colnames(res$fit_stats)))
})

test_that("MetacellsByGroupsSuperCell mcRigor QC populates dubious verdicts", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    skip_if_not(requireNamespace("mcRigor", quietly = TRUE))
    obj <- .sc_test_fixture()
    out <- MetacellsByGroupsSuperCell(
        obj, group.by = "group", min_cells = 50, cells_per_metacell = 10,
        seed = 5, verbose = FALSE
    )
    res <- GetMetacellSuperCell(out)
    expect_true(any(!is.na(res$qc$dubious)))
    expect_true(all(c("mcDiv", "mcDiv_threshold", "zerorate_mc") %in%
                        colnames(res$qc)))
    # QC columns joined onto the metacell metadata
    expect_true("dubious" %in% colnames(res$metacell_object@meta.data))
})

test_that("MetacellsByGroupsSuperCell target_metacells bisects gamma", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    obj <- .sc_test_fixture()
    out <- MetacellsByGroupsSuperCell(
        obj, group.by = "group", min_cells = 50, target_metacells = 10,
        seed = 5, verbose = FALSE
    )
    res <- GetMetacellSuperCell(out)
    n_mc <- ncol(res$metacell_object)
    expect_lte(abs(n_mc - 20), 0.10 * 20)   # 2 groups x ~10
})

test_that("MetacellsByGroupsSuperCell skips small groups and errors when all skip", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    obj <- .sc_test_fixture()
    # mixed: shrink group A below min_cells -> warning, group B still processed
    cells_A <- colnames(obj)[obj$group == "A"]
    cells_B <- colnames(obj)[obj$group == "B"]
    obj_mix <- subset(obj, cells = c(cells_A[1:10], cells_B))
    expect_warning(
        out <- MetacellsByGroupsSuperCell(
            obj_mix, group.by = "group", min_cells = 50,
            cells_per_metacell = 10, seed = 5, verbose = FALSE
        ),
        "skipped"
    )
    expect_equal(nrow(out@misc$SuperCell$skipped_groups), 1L)
    expect_equal(out@misc$SuperCell$skipped_groups$reason, "below min_cells")
    # all groups below min_cells -> hard error
    expect_error(
        MetacellsByGroupsSuperCell(
            obj, group.by = "group", min_cells = 1000, verbose = FALSE
        ),
        "no groups passed min_cells"
    )
})

test_that("OptimizeGammaSuperCell picks a gamma and returns the score table", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    skip_if_not(requireNamespace("mcRigor", quietly = TRUE))
    obj <- .sc_test_fixture()
    res <- OptimizeGammaSuperCell(
        obj, group.by = "group", gammas = c(5, 10, 20), min_cells = 50,
        seed = 5, verbose = FALSE
    )
    expect_equal(length(res$optimal_gamma), 2L)   # both groups
    expect_true(all(res$optimal_gamma %in% c(5, 10, 20)))
    expect_equal(nrow(res$scores), 6L)            # 2 groups x 3 gammas
    expect_true(any(!is.na(res$scores$score)))
    # membership vectors for the chosen granularity, named by cell
    expect_equal(length(res$membership[[1]]), 100L)
    expect_true(all(grepl("^group_", res$membership[[1]])))
})

test_that("mcRigor_DETECT accepts MC-string labels at current HEAD (regression)", {
    skip_if_not(requireNamespace("Seurat", quietly = TRUE))
    skip_if_not(requireNamespace("SuperCell", quietly = TRUE))
    skip_if_not(requireNamespace("mcRigor", quietly = TRUE))
    obj <- .sc_test_fixture()
    counts_g <- .get_counts_matrix(obj, "RNA", "counts", "counts", colnames(obj))
    sc <- SuperCell::SCimplify_for_Seurat_v5(
        obj, gamma = 10, k.knn = 30, assay = "RNA", reduction = list("pca"),
        dims = list(1:5), return.seurat = FALSE, verbose = FALSE
    )
    uid <- sort(unique(sc$membership))
    mem_lab <- paste0("grp_", match(sc$membership, uid))
    names(mem_lab) <- names(sc$membership)
    qc <- cellprograms:::.mcrigor_detect_group(
        counts_g = counts_g, membership = mem_lab, mode = "sum",
        mcrigor_args = list(Nrep = 10), group_id = "A"
    )
    # metacell labels are remapped to the ORIGINAL names, not MC codes
    expect_true(all(qc$metacell %in% mem_lab))
    expect_true(any(!is.na(qc$dubious)))
})
