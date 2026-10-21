# Stage 2a: SuperCell memberships (pure R).
# Usage: Rscript benchmark_memberships_supercell.R <fixture_dir> <out_dir> <scale> <seed>
# Gamma searched per target N to land within +/-10% of the target metacell count.

suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
})

args <- commandArgs(trailingOnly = TRUE)
fx_dir <- args[1]; out_dir <- args[2]; scale <- args[3]; seed <- as.integer(args[4])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(fx_dir, paste0("counts_", scale, "_s", seed, ".rds")))
targets <- c(50, 100, 200, 400)

# Standard preprocessing: normalize -> variable features -> PCA (SuperCell works
# on the cell-to-cell structure; v5 helper handles the Seurat side)
obj <- CreateSeuratObject(counts = counts, assay = "RNA", min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

gamma_for_target <- function(obj, target, tol = 0.10, max_iter = 14) {
  # SCimplify gamma ~= cells per metacell (N_mc = floor(n_cells/gamma));
  # bisection on log-gamma: higher gamma -> fewer metacells
  n_cells <- ncol(obj)
  lo <- max(2, n_cells / (target * 4)); hi <- n_cells / max(1, target / 4)
  best <- NULL
  for (it in seq_len(max_iter)) {
    g <- sqrt(lo * hi)
    sc <- try(SCimplify_for_Seurat_v5(obj, gamma = g, return.seurat = FALSE,
                                      verbose = FALSE), silent = TRUE)
    if (inherits(sc, "try-error") || is.null(sc$membership)) { hi <- g; next }
    n_mc <- length(unique(sc$membership))
    best <- list(gamma = g, n = n_mc, mem = sc$membership)
    if (abs(n_mc - target) <= tol * target) break
    if (n_mc > target) lo <- g else hi <- g
  }
  best
}

for (N in targets) {
  t0 <- Sys.time()
  hit <- gamma_for_target(obj, N)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  write.csv(data.frame(method = "supercell", scale = scale, seed = seed, N = N,
                       minutes = mins),
            file.path(out_dir, sprintf("timing_supercell_%s_s%d_N%d.csv",
                                       scale, seed, N)), row.names = FALSE)
  if (is.null(hit) || is.null(hit$mem)) {
    cat(sprintf("SUPERCELL FAIL %s s%d N=%d: no membership\n", scale, seed, N))
    next
  }
  mem <- hit$mem
  names(mem) <- colnames(obj)
  out <- data.frame(cell_id = names(mem), membership = as.character(mem),
                    stringsAsFactors = FALSE)
  write.csv(out, file.path(out_dir, sprintf("membership_supercell_%s_s%d_N%d.csv",
                                            scale, seed, N)), row.names = FALSE)
  cat(sprintf("SUPERCELL OK %s s%d N=%d: gamma=%.3f got=%d\n",
              scale, seed, N, hit$gamma, hit$n))
}
cat("SUPERCELL STAGE DONE\n")
