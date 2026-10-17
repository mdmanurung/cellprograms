# SENSITIVITY (non-registered): genenmf mc/sc metaprograms at min.confidence=0.2.
# The registered recipe (defaults, min.confidence=0.5) yields 0 MPs on these
# granularities (max consensus confidence 0.03-0.29); this variant quantifies
# what the MP structure looks like if the package confidence bar is lowered.
# Outputs are written to *_conf02.csv and are NOT the registered arm results.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(GeneNMF))

root <- "/mnt/shared-workspace/shared/audit/F_granularity"
MINCONF <- 0.2
for (arm in c("metacell", "singlecell")) {
  allW <- read.csv(sprintf("%s/F_genenmf_%s_allW.csv", root, arm), row.names = 1, check.names = FALSE)
  allH <- read.csv(sprintf("%s/F_genenmf_%s_allH.csv", root, arm), row.names = 1, check.names = FALSE)
  nmf.res <- list()
  for (k in 2:15) {
    Wk <- allW[, grepl(sprintf("^k%d_pattern", k), colnames(allW)), drop = FALSE]
    Hk <- allH[grepl(sprintf("^k%d_pattern", k), rownames(allH)), , drop = FALSE]
    nmf.res[[sprintf("k%d", k)]] <- list(w = as.matrix(Wk), h = as.matrix(Hk))
  }
  mp <- GeneNMF::getMetaPrograms(nmf.res, nMP = 10, min.confidence = MINCONF)
  n_mp <- length(mp$markers.consensus)
  cat(sprintf("[%s] min.confidence=%.2f -> n_mp=%d\n", arm, MINCONF, n_mp))
  write.csv(data.frame(program = names(mp$markers.consensus),
                       genes = vapply(mp$markers.consensus, function(g) paste(g, collapse = ";"), ""),
                       stringsAsFactors = FALSE),
            sprintf("%s/F_genenmf_%s_mp_genes_conf02.csv", root, arm), row.names = FALSE)
  mp$metrics |> write.csv(sprintf("%s/F_genenmf_%s_mp_metrics_conf02.csv", root, arm))
  mp$composition |> write.csv(sprintf("%s/F_genenmf_%s_mp_composition_conf02.csv", root, arm), na = "")
}
cat("sensitivity done\n")
