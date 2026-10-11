## real_tf_triage.R -- annotate_cmp_tf() (CollecTRI) + cmp_triage() on the cap-30 Stephenson and COMBAT CMPs.
## Rscript real_tf_triage.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
pkgload::load_all(".", quiet = TRUE)   # run from the package root
R <- "/mnt/shared-workspace/cmp/real"; O <- "/mnt/results/tusoskill/cmp_v2/interpret"; dir.create(O, FALSE, TRUE)
put <- function(x, f) { tf <- tempfile(fileext = ".csv"); fwrite(x, tf); unlink(file.path(O, f)); system2("cp", c(tf, file.path(O, f))) }

raw <- fread("/mnt/shared-workspace/cmp/priors/collectri_raw.tsv")
net <- prepare_collectri(as.data.frame(raw))
cat(sprintf("CollecTRI: %d TF-target pairs, %d TFs, %d targets\n", nrow(net), uniqueN(net$source), uniqueN(net$target)))
put(net, "collectri_net.csv")

runs <- list(stephenson = list(n = "stephenson_pb_sex_K30", rep = "institute"),
             combat = list(n = "combat_pb_sex_K30", rep = NULL))   # COMBAT: St_Georges has 14 donors (< 15), no replication level
res <- list()
for (co in names(runs)) {
  fit <- readRDS(file.path(R, paste0(runs[[co]]$n, "_fit.rds")))
  cm <- readRDS(file.path(R, paste0(runs[[co]]$n, "_cmp.rds")))$flag
  t0 <- Sys.time()
  tf <- annotate_cmp_tf(fit, cm, net, n_perm = 1000, seed = 1)
  cat(co, format(Sys.time() - t0), "\n"); print(tf)
  tr <- cmp_triage(fit, cm, context = "severity", covariates = c("institute", "sex"),
                   replicate_by = runs[[co]]$rep, tf = tf)
  put(tf$members, sprintf("%s_K30_tf_members.csv", co))
  put(tf$summary, sprintf("%s_K30_tf_summary.csv", co))
  put(tr, sprintf("%s_K30_triage.csv", co))
  res[[co]] <- list(tf = tf, tr = tr)
  print(as.data.table(tr)[, .(cmp, n_members, label, flags, context_p = signif(context_p, 2), min_concordance,
                              tf_shared = substr(tf_shared, 1, 60))])
}

## CMP-level footprint vectors (Stouffer z of member t over all tested TFs) and cross-cohort footprint agreement
zvec <- function(mem) {
  m <- as.data.table(mem)
  m[, z := -sign(t) * qnorm(pt(-abs(t), 1998, log.p = TRUE), log.p = TRUE)]
  dcast(m[, .(z = sum(z) / sqrt(.N)), by = .(cmp, tf)], tf ~ cmp, value.var = "z")
}
ZS <- zvec(res$stephenson$tf$members); ZC <- zvec(res$combat$tf$members)
common <- intersect(ZS$tf, ZC$tf)
A <- as.matrix(ZS[match(common, tf), -1]); B <- as.matrix(ZC[match(common, tf), -1])
C <- cor(A, B, use = "pairwise.complete.obs")
mt <- fread("/mnt/results/tusoskill/cmp_v1/real/cross_cohort_best_match_K30.csv")
print(head(mt))
put(data.table(steph_cmp = rownames(C)[row(C)], combat_cmp = colnames(C)[col(C)], footprint_r = as.vector(C)),
    "cross_cohort_tf_footprint_cor.csv")

## figure: shared TFs per CMP (both cohorts), CMPs not labelled ambient/vdj
fig <- rbindlist(lapply(names(res), function(co) {
  s <- as.data.table(res[[co]]$tf$summary)[class == "shared"]
  lab <- as.data.table(res[[co]]$tr)[, .(cmp, label)]
  s <- merge(s, lab, by = "cmp")
  s[, cohort := co][, cmp_lab := sprintf("%s %d (%s)", ifelse(co == "stephenson", "S", "C"), cmp, label)]
}))
if (nrow(fig)) {
  top_tf <- unique(fig[order(-abs(z))][, head(.SD, 3), by = .(cohort, cmp)]$tf)   # top 3 shared TFs per CMP
  f2 <- fig[tf %in% top_tf]
  f2[, n_sig := pmax(n_up, n_down)]
  f2[, cmp_lab := factor(cmp_lab, levels = unique(f2[order(cohort, cmp), cmp_lab]))]
  ord <- f2[, .(m = max(abs(z))), by = tf][order(m), tf]
  f2[, tf := factor(tf, levels = ord)]
  theme_set(theme_bw(base_family = "Liberation Sans"))
  p <- ggplot(f2, aes(cmp_lab, tf, colour = z, size = n_sig)) + geom_point() +
    scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", name = "Stouffer z") +
    scale_size_continuous(range = c(1.5, 5), name = "Members called") +
    facet_grid(. ~ cohort, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = "Shared TF footprints of cap-30 CMPs (CollecTRI)") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
  tf_png <- tempfile(fileext = ".png")
  ggsave(tf_png, p, width = 12, height = 7, dpi = 150)
  unlink(file.path(O, "fig_shared_tf_footprints_K30.png")); system2("cp", c(tf_png, file.path(O, "fig_shared_tf_footprints_K30.png")))
}
saveRDS(res, "/mnt/shared-workspace/cmp/real/tf_triage_K30.rds")
