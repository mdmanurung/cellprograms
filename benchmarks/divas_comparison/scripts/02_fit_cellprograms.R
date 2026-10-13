# 02_fit_cellprograms.R -- fit cellprograms on the 6 DIVAS COVID-19 blocks
#
# DIVAS-parity recipe (mirrors run_combined_6block_divas.R):
#   - DIVASmain(rowCent = TRUE, colCent = FALSE)  ->  center = TRUE, scale = FALSE
#   - no covariate adjustment (DIVAS ran unadjusted)
#   - all features used (features = "all")
#   cellprograms-specific: point_laplace loading prior, var_type = 1, backfit,
#   nullcheck, max_factors = 30, seed = 1.
#
# Usage: Rscript 02_fit_cellprograms.R <data_dir> <out_dir>

# Repo scripts run against either an installed cellprograms or a source checkout
# (set CELLPROGRAMS_ROOT to the package root for the load_all fallback).
if (requireNamespace("cellprograms", quietly = TRUE)) {
  library(cellprograms)
} else {
  root <- Sys.getenv("CELLPROGRAMS_ROOT", unset = file.path("..", "..", ".."))
  pkgload::load_all(root, quiet = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir  <- if (length(args) >= 2) args[2] else "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BLOCKS <- c(CD4_T = "CD4_T_combined", CD8_T = "CD8_T_combined",
            CD14_Monocyte = "CD14_Monocyte", NK = "NK",
            proteomics = "proteomics", metabolomics = "metabolomics")
GENE_BLOCKS <- c("CD4_T", "CD8_T", "CD14_Monocyte", "NK")

read_block <- function(path) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  stopifnot(is.numeric(mat[1, 1]))
  t(mat)                             # samples x features
}

cat("== loading blocks\n")
pseudobulk <- lapply(BLOCKS, function(b) {
  m <- read_block(file.path(data_dir, "input", paste0("Combined_", b, ".csv")))
  cat("  ", b, ":", paste(dim(m), collapse = " x "), "\n")
  m
})

common <- Reduce(intersect, lapply(pseudobulk, rownames))
cat("common samples:", length(common), "\n")
stopifnot(length(common) == 114L)
pseudobulk <- lapply(pseudobulk, function(m) m[common, , drop = FALSE])

# The 4 gene blocks must share a gene universe for program_sharing (mashr).
gene_sets <- lapply(GENE_BLOCKS, function(b) colnames(pseudobulk[[b]]))
common_genes <- Reduce(intersect, gene_sets)
cat("gene universe: ", length(common_genes),
    " (block sizes: ", paste(vapply(gene_sets, length, 0), collapse = ","), ")\n", sep = "")
for (b in GENE_BLOCKS) pseudobulk[[b]] <- pseudobulk[[b]][, common_genes, drop = FALSE]

meta <- read.csv(file.path(data_dir, "input", "final_sample_info.csv"),
                 stringsAsFactors = FALSE)
meta <- meta[match(common, meta$sample_id), ]
stopifnot(!anyNA(meta$sample_id))
sample_metadata <- data.frame(
  observation_id = meta$sample_id,
  patient   = meta$clinical_id,
  timepoint = meta$timepoint,
  wos       = meta$wos,
  severity  = meta$severity,
  sex       = meta$sex,
  age       = meta$age,
  stringsAsFactors = FALSE
)
cat("severity table:\n"); print(table(sample_metadata$severity))
cat("timepoint table:\n"); print(table(sample_metadata$timepoint))

x <- as_cell_program_data(pseudobulk, sample_metadata = sample_metadata)

cat("== fitting 6 blocks (max_factors=30, point_laplace, backfit, nullcheck, seed=1)\n")
t0 <- Sys.time()
fit <- fit_celltype_programs(
  x,
  loading_prior = "point_laplace",
  max_factors   = 30,
  center        = TRUE,
  scale         = FALSE,
  features      = "all",
  var_type      = 1,
  backfit       = TRUE,
  nullcheck     = TRUE,
  seed          = 1
)
fit <- canonicalize_programs(fit)
cat("total fit time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

k_per_block <- vapply(fit$cell_types, function(ct) ncol(fit$scores[[ct]]), 0)
cat("programs per block:\n"); print(k_per_block)
stopifnot(sum(k_per_block) > 0)

saveRDS(fit, file.path(out_dir, "divas6_fit.rds"))
cat("saved:", file.path(out_dir, "divas6_fit.rds"), "\n")
cat("DONE\n")
