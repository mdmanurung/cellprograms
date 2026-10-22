## quick checks of sim_coord.R (v2: priv_corr, R20small, each-alone)
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/test_sim_coord.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
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

## each-alone families
d <- gen_coord("act_only", amp = 0.5, seed = 11, breadth = "R50")
chk(length(d$truth$programs) == 1 && length(d$truth$programs$R50$members) == 5, "each-alone R50 plants exactly one 5-ct program")
d <- gen_coord("act_only", amp = 0.5, seed = 11, breadth = "R20")
chk(length(d$truth$programs) == 1 && length(d$truth$programs$R20$members) == 2, "each-alone R20 plants exactly one 2-ct program")
d <- gen_coord("act_only", amp = 0.5, seed = 11, breadth = "R20small")
chk(length(d$truth$programs) == 1 && all(d$truth$programs$R20$members %in% d$truth$small), "R20small members are both small cell types")
d <- gen_coord("act+genes", amp = 0.5, seed = 11, breadth = "R100")
chk(length(d$truth$programs) == 1 && length(d$truth$programs$R100$members) == 10, "each-alone R100")

## priv_corr: no recurring programs; 3 decoy pairs with cor ~ cor_priv, disjoint genes
d <- gen_coord("priv_corr", amp = 0.5, seed = 13, cor_priv = 0.7)
tr <- d$truth
chk(length(tr$programs) == 0, "priv_corr plants no recurring programs")
chk(length(tr$decoy_pairs) == 3, "priv_corr plants 3 decoy pairs")
for (dp in tr$decoy_pairs) {
  a <- dp$members[1]; b <- dp$members[2]
  ids <- intersect(names(dp$z[[1]]), names(dp$z[[2]]))
  chk(abs(cor(dp$z[[1]][ids], dp$z[[2]][ids]) - 0.7) < 0.1, paste(dp$id, "activity cor ~0.7"))
  chk(length(intersect(dp$genes[[1]], dp$genes[[2]])) == 0, paste(dp$id, "genes disjoint"))
  chk(!anyDuplicated(dp$members), paste(dp$id, "distinct cell types"))
}
## decoy genes disjoint from private genes
allg <- unlist(lapply(tr$decoy_pairs, function(dp) unlist(dp$genes)))
chk(!anyDuplicated(allg), "priv_corr gene sets mutually disjoint")

## adjusted_matrices: residualized + centered, no flashier needed
d <- gen_coord("act_only", amp = 0.5, seed = 17, breadth = "R50")
mats <- adjusted_matrices(d$x, "institute")
chk(length(mats) == 10 && all(vapply(mats, is.matrix, TRUE)), "adjusted_matrices returns 10 matrices")
chk(all(vapply(mats, function(m) max(abs(colMeans(m))) < 1e-8, TRUE)), "adjusted matrices gene-centered")
## nuisance regression actually removes institute signal
b <- lm(mats$ct01[, 1] ~ d$x$sample_metadata$institute[match(rownames(mats$ct01), d$x$sample_metadata$observation_id)])
chk(sum(coef(b)[-1]^2, na.rm = TRUE) < 1e-8, "institute effect residualized out of gene 1")
cat("all checks passed\n")
