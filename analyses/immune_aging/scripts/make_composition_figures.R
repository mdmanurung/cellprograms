# Composition-level aging figures.
#   fig_aging_composition.png: heatmap of relative proportion change per decade
#     (rows = key cell types, cols = cohorts)
#   fig_aging_tem_ratio.png: GZMK+/GZMB+ Tem log-ratio — OneK1K age scatter +
#     SoundLife CMV boxplot
#
# Usage: Rscript make_composition_figures.R <tables_dir> <out_dir>

suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(gridExtra)
})

args <- commandArgs(TRUE)
tab_dir <- if (length(args) >= 1) args[1] else "/workspace/aging_pbmc/composition"
out_dir <- if (length(args) >= 2) args[2] else "."

traj <- read.csv(file.path(tab_dir, "composition_age_trajectories.csv"), stringsAsFactors = FALSE)
tem <- read.csv(file.path(tab_dir, "tem_ratio_age_cmv.csv"), stringsAsFactors = FALSE)

cohort_order <- c("ABF300", "OneK1K", "AIDA", "CIMA", "KCL", "UCSF", "SoundLife", "SPAC")
type_order <- c("CD8+ RTE", "CD4+ RTE", "CD8+ T Naive", "CD4+ T Naive",
                "gd T Vd2 GZMK+", "CD8+ Tem GZMK+", "CD8+ Tem GZMB+", "CD8+ Temra")
traj$cohort <- factor(traj$cohort, levels = cohort_order)
traj$cell_type <- factor(traj$cell_type, levels = rev(type_order))

p1 <- ggplot(traj, aes(cohort, cell_type, fill = age_rel_change_per10y)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.0f%%", 100 * age_rel_change_per10y)), size = 2.6) +
  scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                       midpoint = 0, name = "Relative change\nper 10y") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Cell-type proportion change with age (prop ~ age + sex)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

# Tem ratio panels
k1k <- list(dir = "onek1k", name = "OneK1K"); sl <- list(dir = "soundlife", name = "SoundLife")
mk_scatter <- function(co) {
  cnt <- read.csv(file.path(tab_dir, paste0(co$dir, "_ct_counts.csv")), stringsAsFactors = FALSE)
  meta <- read.csv(file.path("/workspace/aging_pbmc/pseudobulk", co$dir, "sample_metadata.csv"),
                   stringsAsFactors = FALSE)
  gzk <- cnt[cnt$cell_type == "CD8+ Tem GZMK+", c("donor_id", "n_cells")]
  gzb <- cnt[cnt$cell_type == "CD8+ Tem GZMB+", c("donor_id", "n_cells")]
  d <- merge(merge(gzk, gzb, by = "donor_id", suffixes = c("_gzk", "_gzb")),
             meta[, c("donor_id", "age")], by = "donor_id")
  d$age <- as.numeric(d$age)
  d$ratio <- log((d$n_cells_gzk + 0.5) / (d$n_cells_gzb + 0.5))
  d <- d[complete.cases(d[, c("age", "ratio")]), ]
  ggplot(d, aes(age, ratio)) +
    geom_point(size = 0.9, alpha = 0.55, color = "#0279EE") +
    geom_smooth(method = "lm", color = "#FF9400", linewidth = 0.8, se = FALSE) +
    labs(x = "Age (years)", y = "log GZMK+/GZMB+ ratio",
         title = paste0(co$name, ": age trajectory")) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}
mk_box <- function(co) {
  cnt <- read.csv(file.path(tab_dir, paste0(co$dir, "_ct_counts.csv")), stringsAsFactors = FALSE)
  meta <- read.csv(file.path("/workspace/aging_pbmc/pseudobulk", co$dir, "sample_metadata.csv"),
                   stringsAsFactors = FALSE)
  gzk <- cnt[cnt$cell_type == "CD8+ Tem GZMK+", c("donor_id", "n_cells")]
  gzb <- cnt[cnt$cell_type == "CD8+ Tem GZMB+", c("donor_id", "n_cells")]
  d <- merge(merge(gzk, gzb, by = "donor_id", suffixes = c("_gzk", "_gzb")),
             meta[, c("donor_id", "CMV")], by = "donor_id")
  d$ratio <- log((d$n_cells_gzk + 0.5) / (d$n_cells_gzb + 0.5))
  d <- d[complete.cases(d[, c("CMV", "ratio")]) & d$CMV %in% c("Negative", "Positive"), ]
  d$CMV <- factor(d$CMV, levels = c("Negative", "Positive"))
  ggplot(d, aes(CMV, ratio, fill = CMV)) +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.4) +
    scale_fill_manual(values = c("Negative" = "#75A025", "Positive" = "#FD9BED"), guide = "none") +
    labs(x = "CMV serostatus", y = "log GZMK+/GZMB+ ratio",
         title = paste0(co$name, ": CMV effect")) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

p2 <- mk_scatter(k1k); p3 <- mk_box(sl)
p23 <- arrangeGrob(p2, p3, nrow = 1, widths = c(1.4, 1))

ggsave(file.path(out_dir, "fig_aging_composition.png"), p1,
       width = 7.2, height = 4.2, dpi = 200)
ggsave(file.path(out_dir, "fig_aging_tem_ratio.png"), p23,
       width = 8, height = 3.4, dpi = 200)
cat("figures written to", out_dir, "\n")
