# 06_figures_extra.R -- additional concordance visualizations
#
# Requires the tables written by 05_compare.R:
#   score_cor_matrix.csv                    (DIVAS components x cellprograms scores |cor|)
#   concordance_divas_to_cellprograms.csv   (best matches + best_is_cmp flag)
#   clinical_associations.csv               (per-method wos/severity associations)
#   divas_components.csv                    (component annotation)
#
# Figures:
#   fig_concordance_heatmap.png  -- full |cor| matrix, annotated by class/block
#   fig_breadth_concordance.png  -- best |cor| vs DIVAS pattern breadth (n blocks)
#   fig_effect_agreement.png     -- |rho_wos| agreement for matched pairs (45 deg)
#
# Usage: Rscript 06_figures_extra.R <out_dir>

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "results"

cor_mat <- as.matrix(read.csv(file.path(out_dir, "score_cor_matrix.csv"),
                              row.names = 1, check.names = FALSE))
cor_mat[is.na(cor_mat)] <- 0   # NA cor = no concordance evidence
concord <- read.csv(file.path(out_dir, "concordance_divas_to_cellprograms.csv"),
                    stringsAsFactors = FALSE)
a_all   <- read.csv(file.path(out_dir, "clinical_associations.csv"),
                    stringsAsFactors = FALSE)
dv_annot <- read.csv(file.path(out_dir, "divas_components.csv"),
                     stringsAsFactors = FALSE)
dv_annot <- dv_annot[match(rownames(cor_mat), dv_annot$component), ]
stopifnot(identical(rownames(cor_mat), dv_annot$component))

pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00")
names(pal) <- c("CD4", "CD8", "CD14", "NK", "pro", "metab")
theme_set(theme_minimal(base_family = "Liberation Sans"))

BLOCKS <- c("CD4", "CD8", "CD14", "NK", "pro", "metab")

# ================= fig 1: full concordance heatmap =================
# row annotation: DIVAS class + block-membership strip
blocks_of <- function(pattern) strsplit(pattern, "+", fixed = TRUE)[[1]]
row_blocks <- t(vapply(dv_annot$pattern, function(p) {
  b <- blocks_of(p)
  setNames(as.integer(BLOCKS %in% b), BLOCKS)
}, integer(length(BLOCKS))))
colnames(row_blocks) <- BLOCKS

row_ha <- rowAnnotation(
  class = dv_annot$divas_class,
  n_blocks = anno_barplot(dv_annot$n_blocks, border = FALSE,
                          gp = gpar(fill = "grey75"), height = unit(1.2, "cm")),
  blocks = row_blocks,
  col = list(
    class = c(joint = pal[["CD4"]], partial = pal[["CD8"]], individual = "grey60"),
    blocks = colorRamp2(c(0, 1), c("grey93", "grey25"))
  ),
  annotation_name_side = "top"
)

# column annotation: block of each program (CMP columns -> "CMP")
prog_block <- sub("_.*$", "", colnames(cor_mat))
prog_block[!prog_block %in% BLOCKS] <- "CMP"
col_ha <- columnAnnotation(
  block = factor(prog_block, levels = c(BLOCKS, "CMP")),
  col = list(block = c(pal, CMP = "grey40")),
  annotation_name_side = "right"
)

col_fun <- colorRamp2(c(0, 0.3, 0.6, 1), c("grey95", "#F1E9D8", "#E69F00", "#0072B2"))

ht <- Heatmap(cor_mat,
  name = "|Spearman|",
  col = col_fun,
  top_annotation = col_ha,
  left_annotation = row_ha,
  row_split = factor(dv_annot$divas_class, levels = c("joint", "partial", "individual")),
  column_split = factor(prog_block, levels = c(BLOCKS, "CMP")),
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  row_names_gp = gpar(fontsize = 5),
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "|Spearman|")
)
# NB: the base png() device renders blank output on some headless systems
# (grid output silently dropped); ragg::agg_png is the reliable backend.
ht_file <- file.path(out_dir, "fig_concordance_heatmap.png")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(ht_file, width = 11, height = 9, res = 300, units = "in")
  draw(ht, heatmap_legend_side = "right")
  dev.off()
} else {
  pdf(file.path(out_dir, "fig_concordance_heatmap.pdf"), width = 11, height = 9)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}

# ================= fig 2: concordance vs pattern breadth =================
concord$n_blocks <- dv_annot$n_blocks[match(concord$divas_component, dv_annot$component)]
p2 <- ggplot(concord, aes(x = factor(n_blocks), y = best_abs_cor)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, color = "grey25") +
  geom_jitter(aes(color = divas_class, shape = best_is_cmp),
              width = 0.18, size = 1.7, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(joint = pal[["CD4"]], partial = pal[["CD8"]],
                                individual = "grey60"), name = "DIVAS class") +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                     labels = c(`TRUE` = "CMP consensus", `FALSE` = "single program"),
                     name = "Best cellprograms match") +
  labs(x = "DIVAS pattern breadth (number of blocks sharing the component)",
       y = "Best |Spearman cor| to a cellprograms score",
       title = "Concordance vs sharing breadth",
       subtitle = "Crossbar = median; dashed line = |cor| = 0.5 counterpart threshold") +
  theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out_dir, "fig_breadth_concordance.png"), p2,
       width = 7, height = 4.8, dpi = 300)

# ================= fig 6: clinical effect-size agreement =================
# match()-based lookups: name-indexing with NA would return the whole vector
lk <- function(values, keys, names_vec) values[match(keys, names_vec)]
rho_dv <- lk(a_all$rho_wos[a_all$method == "DIVAS"], concord$divas_component,
             a_all$score[a_all$method == "DIVAS"])
rho_cp <- lk(a_all$rho_wos[a_all$method == "cellprograms"], concord$best_cp_program,
             a_all$score[a_all$method == "cellprograms"])
sig_dv <- lk(a_all$p_wos_bh[a_all$method == "DIVAS"] < 0.05, concord$divas_component,
             a_all$score[a_all$method == "DIVAS"])
sig_cp <- lk(a_all$p_wos_bh[a_all$method == "cellprograms"] < 0.05, concord$best_cp_program,
             a_all$score[a_all$method == "cellprograms"])

agree <- data.frame(
  divas_component = concord$divas_component,
  divas_class = concord$divas_class,
  cp_program = concord$best_cp_program,
  rho_dv = abs(rho_dv),
  rho_cp = abs(rho_cp),
  either_sig = (sig_dv | sig_cp) %in% TRUE,
  stringsAsFactors = FALSE
)
agree <- agree[complete.cases(agree), ]
r_eff <- cor(agree$rho_dv, agree$rho_cp, method = "spearman")

p6 <- ggplot(agree, aes(x = rho_dv, y = rho_cp, color = divas_class,
                        shape = either_sig)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = c(joint = pal[["CD4"]], partial = pal[["CD8"]],
                                individual = "grey60"), name = "DIVAS class") +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                     labels = c(`TRUE` = "wos BH<0.05 (either)", `FALSE` = "n.s."),
                     name = "Clinical signal") +
  coord_fixed(xlim = c(0, max(agree$rho_dv, agree$rho_cp) * 1.05),
              ylim = c(0, max(agree$rho_dv, agree$rho_cp) * 1.05)) +
  labs(x = "|Spearman rho| with wos: DIVAS component",
       y = "|Spearman rho| with wos: matched cellprograms score",
       title = "Clinical effect-size agreement (matched pairs)",
       subtitle = sprintf("Dashed = 45 deg; Spearman rho of effect sizes = %.2f (n = %d pairs)",
                          r_eff, nrow(agree))) +
  theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out_dir, "fig_effect_agreement.png"), p6,
       width = 6.5, height = 5.2, dpi = 300)

cat("DONE: 3 figures written to", out_dir, "\n")
