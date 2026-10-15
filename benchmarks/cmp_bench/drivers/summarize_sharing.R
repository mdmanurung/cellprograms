## summarize_sharing.R -- empirical sharing patterns from the real K30 program_sharing runs.
## Usage: Rscript summarize_sharing.R [in_dir] [cohorts_csv] [out_dir]
##   in_dir   (default /workspace) holding {co}_K30_sharing.rds + {co}_K30_prog_artifact_flags.csv
##   cohorts  (default stephenson,combat)
##   out_dir  (default /mnt/results)
## Outputs: real_K30_sharing_summary.csv (per program, artifact overlay),
## real_K30_sharing_groups.csv (per matched group), real_K30_sharing_patterns.png.

.libPaths("/workspace/.Rlib")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
font <- "Liberation Sans"

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "/workspace"
cos <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else
  c("stephenson", "combat")
out_dir <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "/mnt/results"

summ_all <- list(); grp_all <- list(); pat_all <- list()
for (co in cos) {
  res <- readRDS(file.path(in_dir, sprintf("%s_K30_sharing.rds", co)))
  fl <- read.csv(file.path(in_dir, sprintf("%s_K30_prog_artifact_flags.csv", co)))
  flagcols <- grep("_flag$", names(fl), value = TRUE)
  fl$artifact_flagged <- apply(fl[, flagcols, drop = FALSE], 1, any)

  s <- res$summary
  s$cohort <- co
  s <- merge(s, fl[, c("cell_type", "program_id", "artifact_flagged")],
             by = c("cell_type", "program_id"), all.x = TRUE)
  summ_all[[co]] <- s

  for (gi in seq_along(res$groups)) {
    g <- res$groups[[gi]]
    r <- res$results[[gi]]
    midx <- match(g$program_id, fl$program_id)
    grp_all[[length(grp_all) + 1L]] <- data.frame(
      cohort = co, group_id = gi, n_conditions = nrow(g),
      cell_types = paste(g$cell_type, collapse = ";"),
      programs = paste(g$program_id, collapse = ";"),
      classification = r$classification,
      frac_top_shared = r$frac_top_shared,
      frac_top_private = r$frac_top_private,
      p_all_shared = r$p_all_shared,
      n_artifact_members = sum(fl$artifact_flagged[midx], na.rm = TRUE),
      stringsAsFactors = FALSE)
    pw <- r$pattern_weights
    pat_all[[length(pat_all) + 1L]] <- data.frame(
      cohort = co, group_id = gi, pattern = names(pw), weight = as.numeric(pw),
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, summ_all)
grp <- do.call(rbind, grp_all)
pat <- do.call(rbind, pat_all)

## ---- pattern classes: mashr names singletons by condition; ED_*/PCA_* are
## data-driven; simple_het_* are pairwise-heterogeneity bases; null/identity
## carry "no structure".
pat$klass <- ifelse(pat$pattern == "equal_effects", "equal_effects",
             ifelse(grepl("^(ED|PCA)", pat$pattern), "data-driven",
             ifelse(grepl("^simple_het", pat$pattern), "simple_het",
             ifelse(pat$pattern %in% c("null", "identity"), "null/identity",
                    "singleton"))))
patk <- aggregate(weight ~ cohort + group_id + klass, pat, sum)

## ---- CSVs
write.csv(summ, file.path(out_dir, "real_K30_sharing_summary.csv"), row.names = FALSE)
write.csv(grp, file.path(out_dir, "real_K30_sharing_groups.csv"), row.names = FALSE)

## ---- figure: A) per-group p_all_shared; B) pattern-class weights
grp$label <- paste0(grp$cohort, " G", grp$group_id, " (", grp$n_conditions, ")")
grp <- grp[order(grp$cohort, grp$n_conditions, grp$p_all_shared), ]
grp$label <- factor(grp$label, levels = grp$label)

## panel A: lfsr-based fractions of the top genes per group (the classification basis).
## p_all_shared is uninformative here: with data-driven covariances the canonical
## equal_effects pattern underflows even for strongly shared groups.
gl <- reshape(grp[, c("label", "classification", "frac_top_shared", "frac_top_private")],
              idvar = c("label", "classification"), times = c("shared (>=2 conditions)", "private (1 condition)"),
              v.name = "frac", timevar = "call",
              varying = c("frac_top_shared", "frac_top_private"), direction = "long")
gl$call <- factor(gl$call, levels = c("shared (>=2 conditions)", "private (1 condition)"))
pA <- ggplot(gl, aes(x = frac, y = label, fill = call)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("shared (>=2 conditions)" = "#0279EE",
                               "private (1 condition)" = "#FF9400")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Fraction of top-50 genes called active by lfsr < 0.05",
       y = NULL, fill = NULL) +
  theme_classic(base_family = font, base_size = 9) +
  theme(legend.position = "top", axis.text.y = element_text(size = 6))

patk$label <- factor(
  paste0(patk$cohort, " G", patk$group_id, " (",
         grp$n_conditions[match(paste0(patk$cohort, patk$group_id),
                                paste0(grp$cohort, grp$group_id))], ")"),
  levels = levels(grp$label))
patk$klass <- factor(patk$klass,
                     levels = c("equal_effects", "data-driven", "singleton",
                                "simple_het", "null/identity"))
pB <- ggplot(patk, aes(x = klass, y = label, fill = weight)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1), name = "weight") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = font, base_size = 9) +
  theme(axis.text.y = element_text(size = 6), legend.position = "top",
        panel.grid.major = element_blank())

fig <- pA / pB + plot_annotation(
  title = "Program sharing on real K30 fits",
  theme = theme(text = element_text(family = font, size = 11)))
ggsave(file.path(out_dir, "real_K30_sharing_patterns.png"), fig,
       width = 11, height = max(6, 0.28 * length(levels(grp$label)) + 2),
       dpi = 220, device = ragg::agg_png)

## ---- console digest
cat("== per-cohort digest ==\n")
for (co in cos) {
  s <- summ[summ$cohort == co, ]
  cat(sprintf("%s: %d programs | matched groups: %d | classifications: %s\n", co,
      nrow(s), length(unique(s$group_id[!is.na(s$group_id)])),
      paste(names(table(s$classification)), table(s$classification), sep = "=", collapse = ", ")))
  cat(sprintf("  artifact-flagged programs: %d/%d | flagged among matched: %d/%d\n",
      sum(s$artifact_flagged), nrow(s),
      sum(s$artifact_flagged[!is.na(s$group_id)], na.rm = TRUE), sum(!is.na(s$group_id))))
}
cat("== multi-condition groups (n >= 3) by p_all_shared ==\n")
gm <- grp[grp$n_conditions >= 3, ]
gm <- gm[order(-gm$p_all_shared), ]
if (nrow(gm)) print(gm[1:min(12, nrow(gm)),
                       c("cohort", "group_id", "n_conditions", "classification",
                         "frac_top_shared", "p_all_shared", "n_artifact_members")],
                    digits = 2)
cat("DONE\n")
