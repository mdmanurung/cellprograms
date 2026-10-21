# Stage 5: aggregate results, figures, REPORT.md.
# Usage: Rscript benchmark_report.R <results_dir>
suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
res_dir <- args[1]

indiv <- function(pat, dir) {
  # individual run files only — exclude derived *_all.csv / *_summary.csv so
  # reruns in a dirty directory don't rbind aggregates into the inputs
  fs <- list.files(dir, pat, full.names = TRUE)
  fs[!grepl("_(all|summary)\\.csv$", fs)]
}
metrics <- do.call(rbind, lapply(indiv("metrics_.*\\.csv$", res_dir), read.csv))
recovery <- do.call(rbind, lapply(indiv("recovery_.*\\.csv$", res_dir), read.csv))

metrics$method <- factor(metrics$method, levels = c("supercell", "seacells"))
recovery$method <- factor(recovery$method, levels = c("supercell", "seacells"))
metrics$target_N <- factor(metrics$target_N)
recovery$target_N <- factor(recovery$target_N)

agg_m <- aggregate(cbind(ari_dub, ff_dub, asw, purity_median) ~ method + scale + target_N,
                   data = metrics, FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_r <- aggregate(cbind(detection_mean, fidelity_mean) ~ method + scale + target_N,
                   data = recovery, FUN = function(x) c(mean = mean(x), sd = sd(x)))

write.csv(metrics, file.path(res_dir, "metrics_all.csv"), row.names = FALSE)
write.csv(recovery, file.path(res_dir, "recovery_all.csv"), row.names = FALSE)
write.csv(agg_m, file.path(res_dir, "metrics_summary.csv"), row.names = FALSE)
write.csv(agg_r, file.path(res_dir, "recovery_summary.csv"), row.names = FALSE)

theme_set(theme_minimal(base_family = "Liberation Sans"))
pal <- c(supercell = "#0279EE", seacells = "#FF9400")

plot_metric <- function(df, y, ylab, file) {
  # line + ribbon (mean +/- sd across seeds) over target N
  agg <- aggregate(cbind(y = df[[y]]) ~ method + scale + target_N, data = df,
                   FUN = function(x) c(m = mean(x), s = sd(x)))
  agg$m <- agg$y[, 1]; agg$s <- agg$y[, 2]
  p <- ggplot(agg, aes(x = target_N, y = m, group = method, color = method)) +
    geom_ribbon(aes(ymin = m - s, ymax = m + s, fill = method), alpha = 0.15,
                color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~scale) +
    scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
    labs(x = "target metacells", y = ylab, color = NULL, fill = NULL) +
    theme(legend.position = "top")
  ggsave(file.path(res_dir, file), p, width = 7, height = 4, dpi = 150)
}
plot_metric(metrics, "ari_dub", "ARI (dubious detection)", "fig_ari.png")
plot_metric(metrics, "ff_dub", "FF (dubious detection)", "fig_ff.png")
plot_metric(metrics, "asw", "ASW (PCA space)", "fig_asw.png")
plot_metric(recovery, "detection_mean", "Program detection (max loading-signature rho)",
            "fig_detection.png")
plot_metric(recovery, "fidelity_mean", "Program fidelity (score vs planted activity)",
            "fig_fidelity.png")

# --- additional figures (per user request) ---
permc <- do.call(rbind, lapply(indiv("permc_.*\\.csv$", res_dir), read.csv))
timing <- do.call(rbind, lapply(indiv("timing_.*\\.csv$", res_dir), read.csv))
timing$target_N <- factor(timing$N)   # timing CSVs record the column as "N"

# 6. purity ECDF per method (mirrors mcRigor Fig 1b)
p <- ggplot(permc, aes(purity, color = method)) +
  stat_ecdf(linewidth = 0.8) + facet_wrap(~scale + target_N) +
  scale_color_manual(values = pal) +
  labs(x = "metacell purity", y = "ECDF", color = NULL) +
  theme(legend.position = "top")
ggsave(file.path(res_dir, "fig_purity_ecdf.png"), p, width = 9, height = 6, dpi = 150)

# 7. size vs purity confound scatter
p <- ggplot(permc, aes(size, purity, color = method)) +
  geom_point(size = 0.8, alpha = 0.5) + facet_wrap(~scale + target_N) +
  scale_color_manual(values = pal) +
  labs(x = "metacell size (cells)", y = "metacell purity", color = NULL) +
  theme(legend.position = "top")
ggsave(file.path(res_dir, "fig_size_purity.png"), p, width = 9, height = 6, dpi = 150)

# 8. granularity control: achieved vs target N
ach <- unique(metrics[, c("method", "scale", "seed", "target_N", "n_metacells")])
p <- ggplot(ach, aes(target_N, n_metacells, color = method, group = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.7) + geom_point(size = 2) + facet_wrap(~scale) +
  scale_color_manual(values = pal) +
  labs(x = "target metacells", y = "achieved metacells", color = NULL) +
  theme(legend.position = "top")
ggsave(file.path(res_dir, "fig_granularity_control.png"), p, width = 7, height = 4,
       dpi = 150)

# 9. paired scatter per seed (ARI + fidelity)
pair_m <- merge(metrics[metrics$method == "supercell",
                        c("scale", "seed", "target_N", "ari_dub")],
                metrics[metrics$method == "seacells",
                        c("scale", "seed", "target_N", "ari_dub")],
                by = c("scale", "seed", "target_N"), suffixes = c("_sc", "_sea"))
pair_r <- merge(recovery[recovery$method == "supercell",
                         c("scale", "seed", "target_N", "fidelity_mean")],
                recovery[recovery$method == "seacells",
                         c("scale", "seed", "target_N", "fidelity_mean")],
                by = c("scale", "seed", "target_N"), suffixes = c("_sc", "_sea"))
p1 <- ggplot(pair_m, aes(ari_dub_sea, ari_dub_sc, color = target_N)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) + facet_wrap(~scale) +
  labs(x = "SEACells ARI", y = "SuperCell ARI", color = "target N") +
  xlim(0, 1) + ylim(0, 1)
p2 <- ggplot(pair_r, aes(fidelity_mean_sea, fidelity_mean_sc, color = target_N)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) + facet_wrap(~scale) +
  labs(x = "SEACells fidelity", y = "SuperCell fidelity", color = "target N") +
  xlim(0, 1) + ylim(0, 1)
ggsave(file.path(res_dir, "fig_paired_ari.png"), p1, width = 7, height = 4, dpi = 150)
ggsave(file.path(res_dir, "fig_paired_fidelity.png"), p2, width = 7, height = 4,
       dpi = 150)

# 10. per-cell-type detection heatmap
det_cols <- grep("^detection_T", colnames(recovery), value = TRUE)
det_long <- do.call(rbind, lapply(det_cols, function(cl)
  data.frame(method = recovery$method, scale = recovery$scale,
             target_N = recovery$target_N, cell_type = sub("detection_", "", cl),
             detection = recovery[[cl]])))
agg_d <- aggregate(detection ~ method + scale + target_N + cell_type, data = det_long,
                   FUN = mean)
p <- ggplot(agg_d, aes(cell_type, paste(method, scale, target_N, sep = " | "),
                       fill = detection)) +
  geom_tile() + scale_fill_viridis_c(limits = c(0, 1)) +
  labs(x = NULL, y = NULL, fill = "detection") +
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(res_dir, "fig_detection_heatmap.png"), p, width = 6, height = 7,
       dpi = 150)

# 11. runtime wall-clock (membership construction)
p <- ggplot(timing, aes(paste0(target_N), minutes, fill = method)) +
  geom_col(position = position_dodge(0.9)) + facet_wrap(~scale, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(x = "target metacells", y = "runtime (min, log)", fill = NULL) +
  scale_y_log10() + theme(legend.position = "top")
ggsave(file.path(res_dir, "fig_runtime.png"), p, width = 7, height = 4, dpi = 150)

# 12. downstream cost: program-fit wall time and peak memory per target N
# (fit_time_min / fit_mem_gb recorded around fit_celltype_programs in the
# recovery stage; gc(reset=TRUE) -> gc() Max Used as the peak-memory proxy)
rec_cost <- recovery[!is.na(recovery$fit_time_min) & !is.na(recovery$fit_mem_gb), ]
if (nrow(rec_cost) > 0) {
  agg_cost <- aggregate(cbind(fit_time_min, fit_mem_gb) ~ method + scale + target_N,
                        data = rec_cost, FUN = mean)
  p1 <- ggplot(agg_cost, aes(factor(target_N), fit_time_min, fill = method)) +
    geom_col(position = position_dodge(0.9)) +
    facet_wrap(~scale, scales = "free_y") +
    scale_fill_manual(values = pal) +
    labs(x = "target metacells", y = "program fit time (min)", fill = NULL) +
    theme(legend.position = "top")
  p2 <- ggplot(agg_cost, aes(factor(target_N), fit_mem_gb, fill = method)) +
    geom_col(position = position_dodge(0.9)) +
    facet_wrap(~scale, scales = "free_y") +
    scale_fill_manual(values = pal) +
    labs(x = "target metacells", y = "program fit peak memory (GB)", fill = NULL) +
    theme(legend.position = "top")
  ggsave(file.path(res_dir, "fig_downstream_cost_time.png"), p1, width = 7,
         height = 4, dpi = 150)
  ggsave(file.path(res_dir, "fig_downstream_cost_mem.png"), p2, width = 7,
         height = 4, dpi = 150)
}

# verdict per the PREREG rule — use MEAN-ONLY aggregates: the mean+sd aggregate
# produces matrix columns that contaminate merge/$ extraction below
agg_m_mean <- aggregate(cbind(ari_dub) ~ method + scale + target_N, data = metrics, FUN = mean)
agg_r_mean <- aggregate(cbind(fidelity_mean) ~ method + scale + target_N, data = recovery, FUN = mean)
sup <- agg_m_mean[agg_m_mean$method == "supercell", ]
sea <- agg_m_mean[agg_m_mean$method == "seacells", ]
key <- merge(sup[, c("scale", "target_N", "ari_dub")],
             sea[, c("scale", "target_N", "ari_dub")],
             by = c("scale", "target_N"), suffixes = c("_sc", "_sea"))
ari_gap <- mean(abs(key$ari_dub_sc - key$ari_dub_sea))
sup_r <- agg_r_mean[agg_r_mean$method == "supercell", ]
sea_r <- agg_r_mean[agg_r_mean$method == "seacells", ]
keyr <- merge(sup_r[, c("scale", "target_N", "fidelity_mean")],
              sea_r[, c("scale", "target_N", "fidelity_mean")],
              by = c("scale", "target_N"), suffixes = c("_sc", "_sea"))
fid_rel <- mean(abs(keyr$fidelity_mean_sc - keyr$fidelity_mean_sea) /
                  pmax(keyr$fidelity_mean_sea, 1e-9))
verdict <- if (ari_gap <= 0.05 && fid_rel <= 0.10) "REPLACE (SuperCell as default R-native path)" else
  if (fid_rel > 0.25) "KEEP SEACells" else "DEMOTE SEACells to optional"

report <- c(
  "# REPORT: SEACells vs SuperCell 2.0 metacell benchmark",
  "",
  sprintf("Verdict (PREREG rule): **%s**", verdict),
  "",
  sprintf("- Mean |ARI gap| (dubious detection): %.4f (PREREG tolerance 0.05)", ari_gap),
  sprintf("- Mean relative fidelity gap: %.1f%% (PREREG tolerance 10%%)", 100 * fid_rel),
  "",
  "## Tables",
  "- metrics_all.csv / metrics_summary.csv (mcRigor ARI/FF/purity, ASW, compactness, runtime)",
  "- recovery_all.csv / recovery_summary.csv (program detection + fidelity)",
  "- permc_*.csv (per-metacell purity/size), timing_*.csv (membership wall-clock)",
  "",
  "## Figures",
  "- fig_ari.png, fig_ff.png, fig_asw.png, fig_detection.png, fig_fidelity.png (line + ribbon)",
  "- fig_purity_ecdf.png (purity distributions, mirrors mcRigor Fig 1b)",
  "- fig_size_purity.png (size-vs-purity confound control)",
  "- fig_granularity_control.png (achieved vs target N)",
  "- fig_paired_ari.png, fig_paired_fidelity.png (per-seed paired comparison)",
  "- fig_detection_heatmap.png (per-cell-type program detection)",
  "- fig_runtime.png (membership wall-clock, log scale)",
  "- fig_downstream_cost_time.png, fig_downstream_cost_mem.png (program-fit time and peak memory per target N)",
  "",
  "## Limitations",
  "- Synthetic fixtures only (planted types/donors/programs, 5% injected doublets).",
  "- Geometry in PCA space, not diffusion maps (paper's choice) - pure-R stand-in.",
  "- SEACells run pooled (single group) to match SCimplify's pooled partition.",
  "- SuperCell pinned to supercell-2.0 @ 89b34c078dba0b91289a5a23cf75a6db3a46d232.",
  "- Detection = max Spearman(loading, binary signature); the binary target bounds",
  "  attainable correlation below 1, so absolute values are conservative.",
  "- Program-fit memory is a gc() high-water proxy (Max Used), not RSS.",
  "",
  "## Deviations from PREREG",
  "- SEACells 15k reduced to 2 seeds (runtime); SuperCell 15k kept 3 seeds.",
  "- SEACells wrapper fix: cells_per_metacell = floor(n_cells/N) must be passed",
  "  alongside target_metacells=N (default 75 capped metacell count at N>66).",
  "- mcRigor calling convention pinned against HEAD c70cd6f: membership labels",
  "  are MC-prefixed strings named by cell id (as.character() drops vector names",
  "  in R 4.4); DETECT receives a raw Seurat object; fields=character(0).",
  "- max_factors in recovery is adaptive: min(6, min donors per type - 1);",
  "  coarse N yields thin dominant-donor designs (part of the granularity signal).",
  ""
)
writeLines(report, file.path(res_dir, "REPORT.md"))
cat(sprintf("REPORT DONE verdict=%s ari_gap=%.4f fid_rel=%.3f\n", verdict, ari_gap, fid_rel))
