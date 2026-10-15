## figs_v2.R -- v1.2 held-out figures (seeds 1051-1100). Rscript figs_v2.R <outdir>
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
out <- commandArgs(TRUE)[1]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
E <- "/mnt/shared-workspace/cmp/evals"
POS <- c("cmp_distinct", "cmp_mixed", "context", "ambient", "lr", "weak")
th <- theme_bw(base_size = 11) + theme(text = element_text(family = "Liberation Sans"), panel.grid.minor = element_blank())
pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
save <- function(p, f, w, h) {
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300)
  system2("cp", c(tmp, file.path(out, f)))
  ts <- tempfile(fileext = ".svg")
  on.exit(unlink(ts), add = TRUE)
  svglite::svglite(ts, width = w, height = h); print(p); dev.off()
  system2("cp", c(ts, file.path(out, sub("\\.png$", ".svg", f))))
}
wcsv <- function(d, f) { tmp <- tempfile(fileext = ".csv"); fwrite(d, tmp); system2("cp", c(tmp, file.path(out, f))) }
boot <- function(v, B = 2000) { set.seed(1); b <- replicate(B, mean(sample(v, replace = TRUE))); c(mean(v), quantile(b, c(.025, .975))) }
cp <- function(k, n) binom.test(k, n)$conf.int[1:2]
brd <- function(sig) fifelse(grepl("R20", sig), "R20", fifelse(grepl("R50", sig), "R50", fifelse(sig == "CTX", "CTX", "R100")))
basef <- function(b) factor(b, c("param", "semi"), c("Parametric", "Semi-synthetic (Stephenson)"))

arms <- list("Cap 30 + node-aggregate (winner)" = c("heldout2_v2/K30", "n_agg_noanc"),
             "Cap 10 + node-aggregate (v1.1 finalist)" = c("heldout2_v2/main", "n_agg_noanc"),
             "Baseline: node correlation" = c("heldout2_v2/main", "b_node_cor"),
             "Baseline: principal angles" = c("heldout2_v2/main", "b_pa_cc"))
L <- lapply(arms, function(a) readRDS(file.path(E, a[1], paste0(a[2], ".rds"))))
sig <- rbindlist(lapply(names(L), function(n) as.data.table(L[[n]]$signals)[mode == "marginal" & scenario %in% POS][, arm := n]))
sig[, br := brd(signal)]
sig[, arm := factor(arm, names(arms))]

## 1. F1 by scenario x arm (dataset-level mean F1, bootstrap CI over datasets)
ds <- sig[, .(f1 = mean(f1)), by = .(arm, base, scenario, seed)]
d1 <- ds[, { b <- boot(f1); .(mean = b[1], lo = b[2], hi = b[3], n = .N) }, by = .(arm, base, scenario)]
d1[, scenario := factor(scenario, POS)]
p1 <- ggplot(d1, aes(scenario, mean, colour = arm)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.7), size = 0.25) +
  facet_wrap(~ basef(base), ncol = 1) + scale_colour_manual(values = pal, name = NULL) + th +
  coord_cartesian(ylim = c(0, 1)) + labs(x = NULL, y = "Membership F1 (mean, 95% CI)") +
  theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
save(p1, "fig1_f1_by_scenario_heldout.png", 9, 7); wcsv(d1, "f1_by_scenario_heldout.csv")

## 2. miss rate (signal F1 == 0) by CMP breadth x arm, Clopper-Pearson CI
d2 <- sig[, { k <- sum(f1 == 0); ci <- cp(k, .N); .(miss = k / .N, lo = ci[1], hi = ci[2], n = .N) }, by = .(arm, base, br)]
d2[, br := factor(br, c("R20", "R50", "R100", "CTX"), c("20% of cell types", "50%", "100%", "Context-restricted"))]
p2 <- ggplot(d2, aes(br, miss, colour = arm)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.6), size = 0.25) +
  facet_wrap(~ basef(base), ncol = 2) + scale_colour_manual(values = pal, name = NULL) + th +
  labs(x = "Planted CMP breadth", y = "Miss rate (signal F1 = 0; 95% CP CI)") +
  theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
save(p2, "fig2_miss_by_breadth_heldout.png", 9, 4.5); wcsv(d2, "miss_by_breadth_heldout.csv")

## 3. upstream recall (oracle: best achievable membership given the fitted programs) by fit cap
d3 <- sig[arm %in% names(arms)[1:2], { b <- boot(oracle); .(mean = b[1], lo = b[2], hi = b[3], n = .N) }, by = .(arm, base, br)]
d3[, fit := fifelse(grepl("Cap 30", arm), "Cap 30", "Cap 10")]
d3[, br := factor(br, c("R20", "R50", "R100", "CTX"), c("20%", "50%", "100%", "Context"))]
p3 <- ggplot(d3, aes(br, mean, colour = fit)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.5), size = 0.3) +
  facet_wrap(~ basef(base), ncol = 2) + scale_colour_manual(values = pal[c(2, 1)], name = "Program cap") + th +
  coord_cartesian(ylim = c(0, 1)) + labs(x = "Planted CMP breadth", y = "Upstream recall (oracle F1; 95% CI)")
save(p3, "fig3_upstream_recall_heldout.png", 8, 4); wcsv(d3, "upstream_recall_heldout.csv")

## 4. null FWER: held-out null (50/base) for all arms; pooled with nullx (100/base) where available
nullx <- list("Cap 30 + node-aggregate (winner)" = "nullx_v2/K30comb/n_agg_noanc.rds",
              "Cap 10 + node-aggregate (v1.1 finalist)" = "nullx_v11/main/n_agg_noanc.rds")
d4 <- rbindlist(lapply(names(L), function(nm) {
  s <- as.data.table(L[[nm]]$summary)[scenario == "null" & mode == "marginal", .(base, any_false, src = "held-out")]
  if (!is.null(nullx[[nm]]) && file.exists(file.path(E, nullx[[nm]]))) {
    x <- as.data.table(readRDS(file.path(E, nullx[[nm]]))$summary)[scenario == "null" & mode == "marginal", .(base, any_false, src = "nullx")]
    s <- rbind(s, x)
  }
  rbind(s[src == "held-out", { ci <- cp(sum(any_false), .N); .(set = "Held-out (seeds 1051-1100)", fwer = mean(any_false), lo = ci[1], hi = ci[2], n = .N) }, by = base],
        if (uniqueN(s$src) > 1) s[, { ci <- cp(sum(any_false), .N); .(set = "Held-out + extra null (seeds 31-130)", fwer = mean(any_false), lo = ci[1], hi = ci[2], n = .N) }, by = base])[, arm := nm]
}))
d4[, arm := factor(arm, rev(names(arms)))]
p4 <- ggplot(d4, aes(fwer, arm, colour = basef(base), shape = set)) +
  geom_vline(xintercept = c(0.05, 0.10), linetype = c(3, 2), colour = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.6), size = 0.3) +
  scale_colour_manual(values = pal[1:2], name = "Base") + scale_shape_manual(values = c(16, 17), name = "Null set") + th +
  labs(x = "Family-wise error rate on null data (95% CP CI)", y = NULL) + theme(legend.position = "bottom", legend.box = "vertical")
save(p4, "fig4_null_fwer_heldout.png", 8, 4.5); wcsv(d4, "null_fwer_heldout.csv")

## 5. paired-difference forest: winner minus comparator, pooled and per base
w <- ds[arm == names(arms)[1], .(base, scenario, seed, fw = f1)]
d5 <- rbindlist(lapply(names(arms)[-1], function(nm) {
  m <- merge(w, ds[arm == nm, .(base, scenario, seed, fc = f1)], by = c("base", "scenario", "seed"))
  rbindlist(list(m[, { b <- boot(fw - fc); .(subset = "Pooled", diff = b[1], lo = b[2], hi = b[3], n = .N) }],
                 m[, { b <- boot(fw - fc); .(diff = b[1], lo = b[2], hi = b[3], n = .N) }, by = base][, .(subset = as.character(basef(base)), diff, lo, hi, n)]))[, comparator := nm]
}))
d5[, comparator := factor(comparator, rev(names(arms)[-1]))]
d5[, subset := factor(subset, c("Pooled", levels(basef("param"))))]
p5 <- ggplot(d5, aes(diff, comparator, colour = subset)) + geom_vline(xintercept = 0, colour = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.6), size = 0.3) +
  scale_colour_manual(values = c("#000000", pal[1:2]), name = NULL) + th +
  labs(x = "Paired difference in membership F1, winner minus comparator (95% bootstrap CI)", y = NULL) +
  theme(legend.position = "bottom")
save(p5, "fig5_paired_forest_heldout.png", 8.5, 4); wcsv(d5, "paired_forest_heldout.csv")
cat("figs_v2 done\n")
