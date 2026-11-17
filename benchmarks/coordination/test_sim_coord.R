## quick checks of sim_coord.R (v2: priv_corr, R20small, each-alone)
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/test_sim_coord.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
source(file.path(root, "benchmarks", "coordination", "load_package.R")); load_package(root)
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
## frozen main-run simulator: existing modes reproduce the digests captured before the
## ig_tcr / donor_deletion modes were added (2026-11-17, seed 7, amp 0.5, G = 3000)
golden <- c("act_only" = "210644eef2b4ce0e074f5a4fb21a3f11", "act+genes" = "d3a8bfda5d428aa7917456e52aced46b",
            "null" = "069ea44ca22caa5b03b23d65c552f9dc", "priv_corr" = "17d7180c38c4607b7d3517307f91610e")
for (m in names(golden)) {
  dd <- gen_coord(m, amp = 0.5, seed = 7)
  chk(identical(digest::digest(dd$x$pseudobulk), unname(golden[m])), paste(m, "unchanged by the new modes"))
  chk(is.null(dd$truth$traps), paste(m, "has no traps element"))
}

## donor_deletion: a gene set absent in a subset of donors, in every cell type (a genotype-like trap)
for (a in c("del_frac", "n_del_genes", "del_depth")) {
  args <- list(del_frac = 0.5, n_del_genes = 5L, del_depth = 3)
  args[[a]] <- NULL
  chk(inherits(try(do.call(gen_coord, c(list("donor_deletion", amp = 0.5, seed = 21), args)), silent = TRUE), "try-error"),
      paste("donor_deletion requires", a))
}
d <- gen_coord("donor_deletion", amp = 0.5, seed = 21, del_frac = 0.5, n_del_genes = 5L, del_depth = 3)
tr <- d$truth; tp <- tr$traps$donor_deletion
chk(length(tr$programs) == 0, "donor_deletion plants no recurring programs")
chk(length(tp$genes) == 5L && length(tp$donors) == 50L, "donor_deletion: 5 genes, half of 100 donors")
ct <- "ct01"; Y <- d$x$pseudobulk[[ct]]; del <- intersect(tp$donors, rownames(Y)); keep <- setdiff(rownames(Y), tp$donors)
gn <- colnames(Y)[tp$genes]
chk(mean(Y[del, gn]) < mean(Y[keep, gn]) - 2, "donor_deletion: deleted donors sit about del_depth below the others")
chk(!any(tp$genes %in% unlist(lapply(tr$privates, function(x) unlist(lapply(x, `[[`, "genes"))))), "donor_deletion genes disjoint from private genes")

## ig_tcr: a clonal gene set, large in a few donors, independently per affected cell type
for (a in c("clone_frac", "clone_mult", "n_clone_cts", "n_clone_genes")) {
  args <- list(clone_frac = 0.1, clone_mult = 3, n_clone_cts = 2L, n_clone_genes = 10L)
  args[[a]] <- NULL
  chk(inherits(try(do.call(gen_coord, c(list("ig_tcr", amp = 0.5, seed = 22), args)), silent = TRUE), "try-error"),
      paste("ig_tcr requires", a))
}
d <- gen_coord("ig_tcr", amp = 0.5, seed = 22, clone_frac = 0.1, clone_mult = 3, n_clone_cts = 2L, n_clone_genes = 10L)
tr <- d$truth; ti <- tr$traps$ig_tcr
chk(length(tr$programs) == 0, "ig_tcr plants no recurring programs")
chk(length(ti) == 2L, "ig_tcr affects 2 cell types")
chk(all(vapply(ti, function(x) length(x$genes) == 10L, TRUE)), "ig_tcr: 10 genes per affected cell type")
chk(all(vapply(names(ti), function(ct) abs(length(ti[[ct]]$donors) - 0.1 * length(tr$donors[[ct]])) <= 1, TRUE)), "ig_tcr: about 10% clonal donors per cell type")
chk(!identical(ti[[1]]$donors, ti[[2]]$donors), "ig_tcr: clonal donors drawn independently per cell type")
for (ct in names(ti)) {
  Y <- d$x$pseudobulk[[ct]]; cl <- ti[[ct]]$donors; ot <- setdiff(rownames(Y), cl); gc <- colnames(Y)[ti[[ct]]$genes]
  chk(mean(Y[cl, gc]^2) > 1.3 * mean(Y[ot, gc]^2), paste("ig_tcr: clonal donors carry larger expression in the clonal genes,", ct))
  chk(!any(ti[[ct]]$genes %in% unlist(lapply(tr$privates[[ct]], `[[`, "genes"))), paste("ig_tcr genes disjoint from private genes,", ct))
}
cat("all checks passed\n")
