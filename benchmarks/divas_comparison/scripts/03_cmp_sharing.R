# 03_cmp_sharing.R -- coordination (all 6 blocks) + program sharing (4 gene blocks)
#
# Primary:
#   - find_coordinated_programs, marginal mode, n_perm = 999, seed = 1
#   - program_sharing on the 4 gene blocks (mashr needs a common feature space;
#     proteomics/metabolomics still participate in coordination), n_perm = 199,
#     z_thresh = 1.5 (matches the package's real-data benchmark settings)
# Secondary:
#   - within-severity coordination (mode = "within", context = "severity")
#
# Usage: Rscript 03_cmp_sharing.R <data_dir> <out_dir> [--sensitivity]
#   --sensitivity additionally refits with sex+age residualized and reruns the
#   marginal coordination (headline-call robustness check).

if (requireNamespace("cellprograms", quietly = TRUE)) {
  library(cellprograms)
} else {
  root <- Sys.getenv("CELLPROGRAMS_ROOT", unset = file.path("..", "..", ".."))
  pkgload::load_all(root, quiet = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
sensitivity <- "--sensitivity" %in% args
args <- args[!args %in% "--sensitivity"]
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir  <- if (length(args) >= 2) args[2] else "results"

fit <- readRDS(file.path(out_dir, "divas6_fit.rds"))
GENE_BLOCKS <- c("CD4_T", "CD8_T", "CD14_Monocyte", "NK")

# -- helper: subset a cell_program_fit to a set of blocks (all per-block fields) --
subset_fit_celltypes <- function(fit, keep) {
  stopifnot(all(keep %in% fit$cell_types))
  for (nm in c("pseudobulk", "fits", "canonicalization", "scores", "loadings")) {
    if (!is.null(fit[[nm]])) fit[[nm]] <- fit[[nm]][keep]
  }
  fit$cell_types <- keep
  fit
}

cat("== marginal coordination (6 blocks, n_perm=999)\n")
t0 <- Sys.time()
cmp <- find_coordinated_programs(fit, n_perm = 999, seed = 1, min_donors = 20)
cat("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(cmp, file.path(out_dir, "cmp_marginal.rds"))

cat("== within-severity coordination (6 blocks, n_perm=999)\n")
cmp_within <- find_coordinated_programs(fit, mode = "within", context = "severity",
                                        n_perm = 999, seed = 1, min_donors = 20)
saveRDS(cmp_within, file.path(out_dir, "cmp_within_severity.rds"))

cat("== program sharing (4 gene blocks, n_perm=199, z_thresh=1.5)\n")
t0 <- Sys.time()
gene_fit <- subset_fit_celltypes(fit, GENE_BLOCKS)
sharing <- program_sharing(gene_fit, n_perm = 199, seed = 1, z_thresh = 1.5)
cat("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(sharing, file.path(out_dir, "sharing_geneblocks.rds"))

if (sensitivity) {
  cat("== sensitivity: sex+age residualized refit + marginal coordination\n")
  x <- fit  # cell_program_fit inherits cell_program_data
  fit_resid <- fit_celltype_programs(
    x,
    loading_prior = "point_laplace", max_factors = 30,
    center = TRUE, scale = FALSE, features = "all",
    var_type = 1, backfit = TRUE, nullcheck = TRUE,
    covariates = c("sex", "age"), covariate_mode = "residualize",
    seed = 1
  )
  fit_resid <- canonicalize_programs(fit_resid)
  saveRDS(fit_resid, file.path(out_dir, "divas6_fit_resid.rds"))
  cmp_resid <- find_coordinated_programs(fit_resid, n_perm = 999, seed = 1, min_donors = 20)
  saveRDS(cmp_resid, file.path(out_dir, "cmp_marginal_resid.rds"))
}

cat("DONE\n")
