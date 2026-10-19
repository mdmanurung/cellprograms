# I/O helpers: Seurat -> disk -> Python -> disk -> R round trip.
#
# Follows the export pattern of the hdWGCNA "Alternative metacell algorithms"
# vignette (smorabit.github.io/hdWGCNA/articles/other_metacells.html):
# MatrixMarket + CSV round trips, avoiding h5ad <-> Seurat converters.

.get_counts_matrix <- function(seurat_obj, assay, slot, layer, cells) {
    # Seurat v5 uses layers; v4 uses slots. Support both.
    counts <- tryCatch(
        {
            lays <- SeuratObject::Layers(seurat_obj[[assay]])
            if (length(lays) > 0) {
                lay <- layer
                if (!(lay %in% lays)) {
                    lay <- grep("counts", lays, value = TRUE)[1]
                    if (is.na(lay)) {
                        cli::cli_abort("Assay {.val {assay}} has no layer {.val {layer}} and no counts-like fallback layer; available layers: {.val {lays}}.")
                    }
                    cli::cli_warn("Layer {.val {layer}} not found in assay {.val {assay}}; falling back to layer {.val {lay}}.")
                }
                SeuratObject::LayerData(seurat_obj[[assay]], layer = lay)
            } else {
                NULL
            }
        },
        error = function(e) NULL
    )
    if (is.null(counts)) {
        counts <- SeuratObject::GetAssayData(seurat_obj[[assay]], slot = slot)
    }
    counts[, cells, drop = FALSE]
}

.get_reduction_embeddings <- function(seurat_obj, reduction, dims, cells) {
    emb <- SeuratObject::Embeddings(seurat_obj, reduction)
    if (!is.null(dims)) {
        if (is.character(dims)) {
            missing_d <- setdiff(dims, colnames(emb))
            if (length(missing_d)) {
                cli::cli_abort("{.arg dims} {.val {missing_d}} not found among the {.val {reduction}} reduction's columns ({.val {colnames(emb)}}).")
            }
        } else if (!is.numeric(dims) || any(!is.finite(dims)) ||
                   any(dims < 1) || any(dims > ncol(emb))) {
            cli::cli_abort("{.arg dims} must be column names or valid indices (1..{ncol(emb)}) of the {.val {reduction}} embedding.")
        }
        emb <- emb[, dims, drop = FALSE]
    }
    emb[cells, , drop = FALSE]
}

.export_group_for_python <- function(counts_g, emb_g, out_dir) {
    # counts_g: genes x cells (dgCMatrix); emb_g: cells x dims (matrix) or NULL
    # Cell names that pandas would re-infer as numbers get mangled on the CSV
    # round trip ("01" -> "1", "1e5" -> "100000") and later KeyError in
    # Python; fail closed here with an actionable message.
    cells <- colnames(counts_g)
    num <- suppressWarnings(as.numeric(cells))
    mangled <- !is.na(num) & as.character(num) != cells
    if (any(mangled)) {
        cli::cli_abort("{sum(mangled)} cell name(s) look numeric and would not survive the CSV round trip (e.g. {.val {cells[which(mangled)][1]}}); rename cells with a non-numeric prefix first.")
    }
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    mtx_path <- file.path(out_dir, "counts_genes_by_cells.mtx")
    # CsparseMatrix: the direct dgCMatrix coercion is deprecated (and errors
    # under Matrix.warnDeprecatedCoerce = 2, which flashier's .onLoad sets).
    Matrix::writeMM(as(counts_g, "CsparseMatrix"), mtx_path)
    utils::write.csv(data.frame(cell = colnames(counts_g)),
                     file.path(out_dir, "cells.csv"), row.names = FALSE)
    utils::write.csv(data.frame(gene = rownames(counts_g)),
                     file.path(out_dir, "genes.csv"), row.names = FALSE)
    pca_path <- NULL
    if (!is.null(emb_g)) {
        pca_path <- file.path(out_dir, "pca.csv")
        utils::write.csv(as.data.frame(emb_g), pca_path)
    }
    list(mtx = mtx_path, cells = file.path(out_dir, "cells.csv"),
         genes = file.path(out_dir, "genes.csv"), pca = pca_path)
}

.read_metacell_outputs <- function(out_dir) {
    # metacell matrix: genes x metacells (as written by the python helper)
    # readMM returns dgTMatrix; coerce via CsparseMatrix (the direct dgCMatrix
    # coercion is deprecated and errors under Matrix.warnDeprecatedCoerce = 2).
    mc_mtx <- as(Matrix::readMM(file.path(out_dir, "metacell_mtx.mtx")), "CsparseMatrix")
    mc_obs <- utils::read.csv(file.path(out_dir, "metacell_obs.csv"),
                              row.names = 1, stringsAsFactors = FALSE)
    soft_path <- file.path(out_dir, "soft_assignments.mtx")
    soft <- if (file.exists(soft_path)) {
        as(Matrix::readMM(soft_path), "CsparseMatrix")  # cells x metacells
    } else {
        NULL
    }
    stats_path <- file.path(out_dir, "fit_stats.csv")
    fit_stats <- if (file.exists(stats_path)) {
        utils::read.csv(stats_path, stringsAsFactors = FALSE)
    } else {
        NULL
    }
    list(mtx = mc_mtx, obs = mc_obs, soft = soft, fit_stats = fit_stats)
}
