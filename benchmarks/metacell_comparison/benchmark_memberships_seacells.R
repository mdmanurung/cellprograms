# Stage 2b: SEACells memberships via the package's basilisk path.
# Usage: Rscript benchmark_memberships_seacells.R <fixture_dir> <out_dir> <scale> <seed>
# MUST run in its own R session (reticulate: one Python per session).
# Pooled partition (single group) to match SuperCell's pooled SCimplify comparison.

suppressPackageStartupMessages({
  library(Seurat)
  library(cellprograms)
})

args <- commandArgs(trailingOnly = TRUE)
fx_dir <- args[1]; out_dir <- args[2]; scale <- args[3]; seed <- as.integer(args[4])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(fx_dir, paste0("counts_", scale, "_s", seed, ".rds")))
targets <- c(50, 100, 200, 400)

obj <- CreateSeuratObject(counts = counts, assay = "RNA", min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
# single pooled group so the partition covers all cells (matches SCimplify)
obj$pooled <- "all"
Idents(obj) <- "pooled"

for (N in targets) {
  t0 <- Sys.time()
  ok <- tryCatch({
    mc <- MetacellsByGroupsSEACells(obj, group.by = "pooled",
                                    target_metacells = N,
                                    cells_per_metacell = max(2, floor(nrow(obj) / N)),
                                    min_cells = 1, verbose = FALSE)
    gm <- GetMetacellSEACells(mc)
    soft <- gm$soft_assignments
    if (is.null(soft)) stop("no soft_assignments in metacell result")
    # hard label = row-wise argmax over the joint soft-assignment matrix
    # (same reconstruction the package uses for mcRigor membership)
    j <- NULL  # (argmax above; no scalar needed)
    mem_vec <- colnames(soft)[max.col(soft, ties.method = "first")]
    names(mem_vec) <- rownames(soft)
    out <- data.frame(cell_id = names(mem_vec), membership = mem_vec,
                      stringsAsFactors = FALSE)
    write.csv(out, file.path(out_dir, sprintf("membership_seacells_%s_s%d_N%d.csv",
                                              scale, seed, N)), row.names = FALSE)
    cat(sprintf("SEACELLS OK %s s%d N=%d: got=%d\n", scale, seed, N,
                length(unique(mem_vec))))
    TRUE
  }, error = function(e) {
    cat(sprintf("SEACELLS FAIL %s s%d N=%d: %s\n", scale, seed, N,
                conditionMessage(e)))
    FALSE
  })
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  write.csv(data.frame(method = "seacells", scale = scale, seed = seed, N = N,
                       minutes = mins),
            file.path(out_dir, sprintf("timing_seacells_%s_s%d_N%d.csv",
                                       scale, seed, N)), row.names = FALSE)
  invisible(ok)
}
cat("SEACELLS STAGE DONE\n")
