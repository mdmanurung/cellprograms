## figs.R -- benchmark figures. Rscript figs.R <outdir>
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
out <- commandArgs(TRUE)[1]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
E <- "/mnt/shared-workspace/cmp/evals"
th <- theme_bw(base_size = 11) + theme(text = element_text(family = "Liberation Sans"), panel.grid.minor = element_blank())
pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#999999", "#F0E442")
save <- function(p, f, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 200); system2("cp", c(tmp, file.path(out, f))) }
rd <- function(d, cn) { f <- file.path(E, d, paste0(cn, ".rds")); if (file.exists(f)) readRDS(f) else NULL }
per_ds <- function(x) { g <- as.data.table(x$signals)[mode == "marginal"]; g[, .(f1 = mean(f1)), by = .(base, scenario, seed)] }
boot <- function(v) { set.seed(1); b <- replicate(2000, mean(sample(v, replace = TRUE))); c(mean(v), quantile(b, c(.025, .975))) }
lab <- c(n_agg_noanc = "Finalist (node-aggregate)", b_node_cor = "Baseline: node correlation", b_pa_cc = "Baseline: principal angles",
         e_rf_maxT_gr_amb3 = "Edge max-T", n_agg_t2 = "Node-aggregate + anchors", n_agg_minp_lt35 = "Node+edge minP",
         b_gene_match = "Gene matching", s_shuffled = "Shuffled control")

## 1. per-scenario F1 (split given by env FIG_SPLIT, default test_v11)
sp <- Sys.getenv("FIG_SPLIT", "test_v11")
cands <- names(lab)
d1 <- rbindlist(lapply(cands, function(cn) { x <- rd(file.path(sp, "main"), cn); if (is.null(x)) return(NULL)
  p <- per_ds(x); p[, { b <- boot(f1); .(mean = b[1], lo = b[2], hi = b[3], n = .N) }, by = .(base, scenario)][, candidate := cn] }))
if (nrow(d1)) {
  d1[, candidate := factor(lab[candidate], lab)]
  d1[, base := factor(base, c("param", "semi"), c("Parametric", "Semi-synthetic (Stephenson)"))]
  p1 <- ggplot(d1, aes(scenario, mean, colour = candidate)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), size = 0.25) +
    facet_wrap(~base, ncol = 1) + scale_colour_manual(values = pal, name = NULL) + th +
    labs(x = NULL, y = "Membership F1 (mean, 95% CI)") + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 3))
  save(p1, sprintf("fig1_f1_by_scenario_%s.png", sp), 9, 7)
  fwrite(d1, file.path(tempdir(), "d1.csv")); system2("cp", c(file.path(tempdir(), "d1.csv"), file.path(out, sprintf("f1_by_scenario_%s.csv", sp))))
}

## 2. null FWER with Clopper-Pearson CI: nullx (seeds 31-130) + test null (1001-1050) where available
d2 <- rbindlist(lapply(cands, function(cn) {
  s <- rbindlist(lapply(c("nullx_v11", sp), function(d) { x <- rd(file.path(d, "main"), cn); if (is.null(x)) NULL else
    as.data.table(x$summary)[scenario == "null" & mode == "marginal", .(base, any_false, src = d)] }))
  if (!nrow(s)) return(NULL)
  s[, { k <- sum(any_false); n <- .N; ci <- binom.test(k, n)$conf.int; .(fwer = k / n, lo = ci[1], hi = ci[2], n = n) }, by = base][, candidate := cn]
}))
if (nrow(d2)) {
  d2[, candidate := factor(lab[candidate], rev(lab))]
  p2 <- ggplot(d2, aes(fwer, candidate, colour = base)) + geom_vline(xintercept = c(0.05, 0.10), linetype = c(3, 2), colour = "grey40") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.5), size = 0.3) +
    scale_colour_manual(values = pal[c(1, 2)], name = "Base") + th + labs(x = "Family-wise error rate on null data (95% CP CI)", y = NULL)
  save(p2, "fig2_null_fwer.png", 7, 4)
  fwrite(d2, file.path(tempdir(), "d2.csv")); system2("cp", c(file.path(tempdir(), "d2.csv"), file.path(out, "null_fwer.csv")))
}

## 3. ablation forest: paired difference vs finalist on dev (same arm), plus fit-arm ablations
fin <- per_ds(rd("dev_v11/main", "n_agg_noanc"))
abl <- list(
  "Add edge anchors (n_agg_t2)" = c("dev_v11/main", "n_agg_t2"),
  "Anchors + link 3.5" = c("dev_v11/main", "n_agg_lt35"),
  "Node+edge minP" = c("dev_v11/main", "n_agg_minp_lt35"),
  "Edge max-T, rank-first" = c("dev_v11/main", "e_rf_maxT_gr_amb3"),
  "Edge max-T, residualize-first" = c("dev_v11/main", "e_rg_maxT_gr_amb3"),
  "Node-aggregate, residualize-first" = c("dev_v11/main", "n_agg_rg"),
  "No ambient filter (anchors)" = c("dev_v11/main", "n_agg_noamb"),
  "Baseline: node correlation" = c("dev_v11/main", "b_node_cor"),
  "Baseline: principal angles" = c("dev_v11/main", "b_pa_cc"),
  "Fit without site adjustment" = c("dev_v11/noadj", "n_agg_noanc"),
  "Fit with K<=20" = c("dev_v11/K20", "n_agg_noanc"),
  "Fit on all genes" = c("dev_v11/featall", "n_agg_noanc"))
d3 <- rbindlist(lapply(names(abl), function(k) { x <- rd(abl[[k]][1], abl[[k]][2]); if (is.null(x)) return(NULL)
  m <- merge(fin, per_ds(x), by = c("base", "scenario", "seed")); b <- boot(m$f1.y - m$f1.x)
  data.table(ablation = k, diff = b[1], lo = b[2], hi = b[3], n = nrow(m)) }))
if (nrow(d3)) {
  d3[, ablation := factor(ablation, rev(names(abl)))]
  p3 <- ggplot(d3, aes(diff, ablation)) + geom_vline(xintercept = 0, colour = "grey40") +
    geom_pointrange(aes(xmin = lo, xmax = hi), colour = pal[1], size = 0.3) + th +
    labs(x = "Paired difference in membership F1 vs finalist (95% bootstrap CI)", y = NULL)
  save(p3, "fig3_ablation_forest.png", 7.5, 4.5)
  fwrite(d3, file.path(tempdir(), "d3.csv")); system2("cp", c(file.path(tempdir(), "d3.csv"), file.path(out, "ablation_paired.csv")))
}
cat("figs done\n")
