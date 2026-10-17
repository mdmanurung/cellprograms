.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(GeneNMF))

root <- "/mnt/shared-workspace/shared/audit/F_granularity"
arm <- "metacell"
allW <- read.csv(sprintf("%s/F_genenmf_%s_allW.csv", root, arm), row.names = 1, check.names = FALSE)
allH <- read.csv(sprintf("%s/F_genenmf_%s_allH.csv", root, arm), row.names = 1, check.names = FALSE)
cat("allW:", dim(allW)[1], "x", dim(allW)[2], " allH:", dim(allH)[1], "x", dim(allH)[2], "\n")

# rebuild minimal nmf.res: list of k2..k15 with $w (genes x k), $h (k x obs)
nmf.res <- list()
for (k in 2:15) {
  wc <- sprintf("^k%d_pattern", k)
  Wk <- allW[, grepl(wc, colnames(allW)), drop = FALSE]
  Hk <- allH[grepl(wc, rownames(allH)), , drop = FALSE]
  nmf.res[[sprintf("k%d", k)]] <- list(w = as.matrix(Wk), h = as.matrix(Hk))
}
wl <- getFromNamespace("weightedLoadings", "GeneNMF")
gmc <- getFromNamespace("get_metaprogram_consensus", "GeneNMF")
gnf <- getFromNamespace("getNMFgenes", "GeneNMF")
jsim <- getFromNamespace("jaccardSimilarity", "GeneNMF")

nmf.wgt <- wl(nmf.res, specificity.weight = 5)
genes <- gnf(nmf.res = nmf.wgt, specificity.weight = NULL, weight.explained = 0.5, max.genes = 200)
cat("patterns with gene sets:", length(genes$nmf.genes), "\n")
J <- jsim(genes)
Jdist <- as.dist(1 - J)
tree <- hclust(Jdist, method = "average")
cl <- cutree(tree, k = 10)
print(table(cl))

# consensus: inspect confidence distribution per cluster
cons <- gmc(nmf.wgt = nmf.wgt, nMP = 10, min.confidence = 0, weight.explained = 0.5,
            max.genes = 200, cl_members = cl)
ng <- vapply(cons, function(m) sum(m >= 0.5), integer(1))
ng0 <- vapply(cons, length, integer(1))
cat("genes per cluster (any conf):", paste(ng0, collapse = ","), "\n")
cat("genes per cluster (conf>=0.5):", paste(ng, collapse = ","), "\n")
cat("max confidence per cluster:  ", paste(vapply(cons, function(m) round(max(c(m, 0)), 2), numeric(1)), collapse = ","), "\n")
