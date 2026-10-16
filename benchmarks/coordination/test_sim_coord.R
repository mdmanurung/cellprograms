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
## artifact decoys: vdj (Ig/TCR-like) and deletion (donor spike)
d <- gen_coord("act+genes", amp = 0.5, seed = 11, artifacts = c("vdj", "deletion"))
tr <- d$truth
chk(setequal(names(tr$programs), c("R100", "R50", "R20", "VDJ", "DEL")), "artifacts planted alongside targets")
chk(identical(tr$programs$VDJ$type, "artifact_vdj") && identical(tr$programs$DEL$type, "artifact_deletion"), "artifact types recorded")
vdj_cols <- grep("^IG[KL]V|^TRBV", colnames(d$x$pseudobulk$ct01), value = TRUE)
chk(length(vdj_cols) == 40, "40 receptor-named genes in the pseudobulk")
chk(length(cellprograms::artifact_genes(vdj_cols, sets = "vdj")) == 40, "package artifact regex flags all 40")
chk(length(tr$programs$VDJ$members) %in% 3:5, "vdj program spans 3-5 cell types")
zv <- tr$programs$VDJ$z; m1 <- names(zv)[1]
zc <- sapply(names(zv)[-1], function(ct) cor(zv[[m1]], zv[[ct]]))
chk(all(zc > 0.999), "vdj z shared across members (composition latent)")
zd <- tr$programs$DEL$z[[1]]
chk(all(tr$programs$DEL$members == COORD_CTS), "deletion acts in every cell type")
chk(setequal(names(zd)[zd == 1], names(zd)[zd == 1]) && sum(zd) %in% 1:2, "deletion touches 1-2 donors")
## targets unaffected: R50 genes still disjoint from the artifact gene sets
art_g <- unlist(lapply(tr$programs[c("VDJ", "DEL")], function(p) unlist(p$genes)))
chk(length(intersect(unlist(tr$programs$R50$genes), art_g)) == 0, "artifact genes disjoint from targets")
## scoring: a VDJ-majority cluster is an artifact call, not a target
source(file.path(here, "methods.R")); source(file.path(here, "scoring.R"))
fitv <- list(cell_types = COORD_CTS, scores = lapply(d$x$pseudobulk, function(Y) NULL))
cl_v <- data.frame(cluster = 1, cell_type = tr$programs$VDJ$members,
                   program_id = paste0(tr$programs$VDJ$members, "_1"), stat = 5, fwer_p = 0.01,
                   stringsAsFactors = FALSE)
## build a fake fit whose first node of each VDJ member is the vdj latent itself
sc <- list()
for (ct in COORD_CTS) {
  Z <- matrix(rnorm(100 * 3), 100, 3, dimnames = list(paste0("d", sprintf("%03d", 1:100)), paste0(ct, "_", 1:3)))
  if (ct %in% tr$programs$VDJ$members) Z[, 1] <- 2 * tr$programs$VDJ$z[[ct]][rownames(Z)]
  sc[[ct]] <- Z
}
fitv <- list(cell_types = COORD_CTS, scores = sc)
ev <- eval_calls_coord(cl_v, fitv, tr, d$strata)
chk(ev$summary$n_artifact == 1 && !isTRUE(ev$summary$any_decoy), "VDJ-majority cluster labelled artifact")
chk(all(ev$f1 == 0), "artifact cluster contributes no F1")

## count-level variant: integer counts, log2 CPM fit input, truth carried over
dc <- gen_coord_counts("act+genes", amp = 0.5, seed = 13, artifacts = c("vdj", "deletion"))
chk(all(unlist(lapply(dc$counts, function(Y) Y >= 0)) ), "counts nonnegative")
chk(all(dc$counts$ct01 == round(dc$counts$ct01)), "counts are integers")
chk(abs(sum(dc$counts$ct01[1, ]) - dc$ncell[["ct01"]][["d001"]]) < 5 * sqrt(dc$ncell[["ct01"]][["d001"]]),
    "donor totals track library sizes (Poisson)")
chk(!identical(dc$x$pseudobulk$ct01, dc$counts$ct01), "fit input is log2 CPM, not counts")
chk(setequal(names(dc$truth$programs), c("R100", "R50", "R20", "VDJ", "DEL")), "count variant carries the truth")
cat("all checks passed\n")
