# Top genes for the strongest-age program of key cell types (OneK1K, largest cohort).
suppressPackageStartupMessages(library(cellprograms))
fit <- readRDS("/mnt/shared-workspace/aging_pbmc/pseudobulk/onek1k/onek1k_full_fit.rds")
assoc <- read.csv("/mnt/shared-workspace/aging_pbmc/pseudobulk/onek1k/onek1k_full_assoc.csv", stringsAsFactors = FALSE)
cts <- c("CD8+ RTE", "CD4+ RTE", "CD8+ Tem GZMK+", "CD8+ Tem GZMB+", "CD8+ T Naive", "CD4+ T Naive")
for (ct in cts) {
  s <- assoc[assoc$cell_type == ct & !is.na(assoc$age_p_fdr), ]
  if (!nrow(s)) { cat(ct, ": no programs tested\n"); next }
  w <- s[which.min(s$age_p_fdr), ]
  prog <- w$program
  k <- sub(".*_", "", prog)  # program index = suffix after last underscore
  f <- fit$fits[[ct]]
  W <- f$flash$F_pm
  ki <- suppressWarnings(as.integer(k))
  if (is.na(ki) || ki > ncol(W)) { cat(ct, ": program index not found\n"); next }
  wv <- W[, ki]
  names(wv) <- rownames(W)
  top <- names(sort(wv, decreasing = TRUE))[1:8]
  bot <- names(sort(wv, decreasing = FALSE))[1:8]
  cat("\n==", ct, "|", prog, "| age beta =", round(w$age_beta_per10y_sd, 3), "SD/10y, FDR =", signif(w$age_p_fdr, 2), "\n")
  cat("  + loadings:", paste(top, collapse = ", "), "\n")
  cat("  - loadings:", paste(bot, collapse = ", "), "\n")
}
