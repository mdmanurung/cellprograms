# Build EBMF preprocessing arms with scran HVG selection per cell type.
# Usage: Rscript preprocess_arms.R <pb_sum_csv_dir> <pb_mean_csv_dir> <outdir>
# Expects per-cell-type CSVs <ct>.csv (samples x genes) in each input dir.
suppressMessages({
  .libPaths(c("/workspace/.Rlib", .libPaths()))
  library(edgeR); library(scran)
})
args <- commandArgs(trailingOnly = TRUE)
sum_dir <- args[1]; mean_dir <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cts <- sub("\\.csv$", "", list.files(sum_dir, pattern = "\\.csv$"))
for (ct in cts) {
  counts <- read.csv(file.path(sum_dir, paste0(ct, ".csv")), row.names = 1, check.names = FALSE)
  # TMM normalization on raw-count sums -> log1p CPM (MOFA-B recipe)
  dge <- edgeR::DGEList(counts = t(counts))  # edgeR wants genes x samples
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  sfs <- dge$samples$lib.size * dge$samples$norm.factors
  tmm_log <- log1p(t(t(counts) / sfs) * 1e6)  # samples x genes

  # arm 1: TMM + scran top-2000 HVGs (selected on the TMM-log matrix)
  keep1 <- scran::getTopHVGs(scran::modelGeneVar(t(tmm_log)), n = 2000)  # genes as rows
  write.csv(tmm_log[, keep1, drop = FALSE],
            file.path(outdir, sprintf("tmm_hvg__%s.csv", ct)), row.names = TRUE)

  # arm 2: log-norm mean pseudobulk + scran top-2000 HVGs
  mean_mat <- read.csv(file.path(mean_dir, paste0(ct, ".csv")), row.names = 1, check.names = FALSE)
  keep2 <- scran::getTopHVGs(scran::modelGeneVar(t(as.matrix(mean_mat))), n = 2000)
  write.csv(mean_mat[, keep2, drop = FALSE],
            file.path(outdir, sprintf("lognorm_hvg__%s.csv", ct)), row.names = TRUE)
  cat("done:", ct, "\n")
}
