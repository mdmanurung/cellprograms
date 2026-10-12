## cross_cohort_v2.R -- as cross_cohort.R, but run names via args: Rscript cross_cohort_v2.R <steph_run> <combat_run> [suffix]
## For CMP i (Stephenson) and j (COMBAT): mean over shared cell types of the Pearson correlation of the members'
## signed loadings (sign = member sign within its CMP) over shared genes. Also best-match table.
suppressPackageStartupMessages(library(data.table))
R <- "/mnt/shared-workspace/cmp"; O <- "/mnt/results/tusoskill/cmp_v1/real"
ld <- function(n) list(fit = readRDS(file.path(R, "real", paste0(n, "_fit.rds"))),
                       cm = as.data.table(readRDS(file.path(R, "real", paste0(n, "_cmp.rds")))$flag$cmps),
                       s = fread(file.path(O, n, "flag_cmp_summary.csv")))
a_ <- commandArgs(TRUE); an <- if (length(a_) >= 2) a_[1] else "stephenson_pb"; bn <- if (length(a_) >= 2) a_[2] else "combat_pb"; sfx <- if (length(a_) >= 3) a_[3] else ""
A <- ld(an); B <- ld(bn)
vec <- function(X, ct, pr, sg) { w <- X$fit$loadings[[ct]][, pr] * sg; w }
M <- matrix(NA_real_, uniqueN(A$cm$cmp), uniqueN(B$cm$cmp), dimnames = list(paste0("S", sort(unique(A$cm$cmp))), paste0("C", sort(unique(B$cm$cmp)))))
NCT <- M
for (i in sort(unique(A$cm$cmp))) for (j in sort(unique(B$cm$cmp))) {
  a <- A$cm[cmp == i]; b <- B$cm[cmp == j]; cts <- intersect(a$cell_type, b$cell_type)
  if (!length(cts)) next
  r <- vapply(cts, function(ct) {
    wa <- vec(A, ct, a[cell_type == ct]$program_id, a[cell_type == ct]$sign)
    wb <- vec(B, ct, b[cell_type == ct]$program_id, b[cell_type == ct]$sign)
    g <- intersect(names(wa), names(wb)); if (length(g) < 50) return(NA_real_)
    stats::cor(wa[g], wb[g]) }, 1)
  ## the CMP sign convention is arbitrary per CMP, so align by the sign of the mean
  M[paste0("S", i), paste0("C", j)] <- abs(mean(r, na.rm = TRUE)); NCT[paste0("S", i), paste0("C", j)] <- sum(!is.na(r))
}
best <- data.table(steph_cmp = rownames(M), combat_best = colnames(M)[apply(M, 1, which.max)], mean_abs_r = round(apply(M, 1, max, na.rm = TRUE), 3),
                   n_shared_ct = NCT[cbind(seq_len(nrow(M)), apply(M, 1, which.max))])
best[, steph_genes := substr(A$s$recurrent_top_genes[match(as.integer(sub("S", "", steph_cmp)), A$s$cmp)], 1, 60)]
best[, combat_genes := substr(B$s$recurrent_top_genes[match(as.integer(sub("C", "", combat_best)), B$s$cmp)], 1, 60)]
print(best)
fwrite(best, "/tmp/cc_best.csv"); fwrite(as.data.table(M, keep.rownames = "steph_cmp"), "/tmp/cc_mat.csv")
system2("cp", c("/tmp/cc_best.csv", file.path(O, paste0("cross_cohort_best_match", sfx, ".csv")))); system2("cp", c("/tmp/cc_mat.csv", file.path(O, paste0("cross_cohort_matrix", sfx, ".csv"))))
