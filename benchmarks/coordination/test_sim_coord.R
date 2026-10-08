## quick checks of sim_coord.R: R_LIBS_USER=/nonexistent Rscript benchmarks/coordination/test_sim_coord.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
source(file.path(here, "sim_coord.R"))
options(ms.G = 3000L)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")
for (mode in c("act+genes", "act+half", "act_only", "genes_only", "null")) {
  d <- gen_coord(mode, amp = 0.5, seed = 7)
  tr <- d$truth
  chk(length(d$x$cell_types) == 10, paste(mode, "10 cell types"))
  chk(all(tr$overlap[tr$small, tr$small][upper.tri(diag(3))] >= 30), paste(mode, "small cell types share >= 30 donors"))
  if (mode == "null") { chk(length(tr$programs) == 0, "null has no recurring programs"); next }
  chk(identical(sapply(tr$programs, function(p) length(p$members)), c(R100 = 10L, R50 = 5L, R20 = 2L)), paste(mode, "breadth 10/5/2"))
  pr <- tr$programs$R50; m <- pr$members
  zc <- sapply(m[-1], function(ct) { a <- intersect(names(pr$z[[m[1]]]), names(pr$z[[ct]])); cor(pr$z[[m[1]]][a], pr$z[[ct]][a]) })
  if (mode == "genes_only") chk(all(abs(zc) < 0.5), "genes_only: z independent across members") else chk(all(zc > 0.999), paste(mode, "z shared across members"))
  g1 <- pr$genes[[m[1]]]; g2 <- pr$genes[[m[2]]]
  ov <- length(intersect(g1, g2)) / length(g1)
  chk(switch(mode, "act+genes" = ov == 1, genes_only = ov == 1, "act+half" = abs(ov - 0.5) < 0.05, act_only = ov == 0), paste(mode, "gene overlap", ov))
  other <- unlist(lapply(tr$programs[c("R100", "R20")], function(p) unlist(p$genes)))
  chk(length(intersect(unlist(pr$genes), other)) == 0, paste(mode, "R50 genes disjoint from R100/R20"))
}
cat("all checks passed\n")
