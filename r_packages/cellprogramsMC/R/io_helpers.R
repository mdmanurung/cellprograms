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
        emb <- emb[, dims, drop = FALSE]
    }
    emb[cells, , drop = FALSE]
}

.export_group_for_python <- function(counts_g, emb_g, out_dir) {
    # counts_g: genes x cells (dgCMatrix); emb_g: cells x dims (matrix) or NULL
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    mtx_path <- file.path(out_dir, "counts_genes_by_cells.mtx")
    Matrix::writeMM(as(counts_g, "dgCMatrix"), mtx_path)
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
    mc_mtx <- as(Matrix::readMM(file.path(out_dir, "metacell_mtx.mtx")), "dgCMatrix")
    mc_obs <- utils::read.csv(file.path(out_dir, "metacell_obs.csv"),
                              row.names = 1, stringsAsFactors = FALSE)
    soft_path <- file.path(out_dir, "soft_assignments.mtx")
    soft <- if (file.exists(soft_path)) {
        as(Matrix::readMM(soft_path), "dgCMatrix")  # cells x metacells
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
