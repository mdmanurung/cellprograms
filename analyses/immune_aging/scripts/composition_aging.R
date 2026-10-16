# Composition-level aging axis: key-type proportions and the GZMK+/GZMB+ Tem ratio.
#
# Per cohort (counts from composition_counts.py, metadata from pseudobulk
# sample_metadata.csv):
#   - key-type proportions: prop ~ age10 + sex [+CMV]; relative change per 10y
#     = beta*1 / mean_prop (beta is per decade because age10 = age/10)
#   - Tem ratio: log((nGZMK+.5)/(nGZMB+.5)) ~ age10 + sex [+CMV]
# CMV term only where CMV is informative (KCL: Neg/Pos/No_info; SoundLife:
# Neg/Pos); elsewhere the column is all-No_info and is dropped.
#
# Usage: Rscript composition_aging.R <composition_dir> <pseudobulk_base> <out_dir>

args <- commandArgs(TRUE)
comp_dir <- if (length(args) >= 1) args[1] else "/workspace/aging_pbmc/composition"
pb_base  <- if (length(args) >= 2) args[2] else "/workspace/aging_pbmc/pseudobulk"
out_dir  <- if (length(args) >= 3) args[3] else comp_dir

cohorts <- list(
  list(name = "ABF300", dir = "abf300"), list(name = "OneK1K", dir = "onek1k"),
  list(name = "AIDA", dir = "aida"), list(name = "CIMA", dir = "cima"),
  list(name = "KCL", dir = "kcl"), list(name = "UCSF", dir = "ucsf"),
  list(name = "SoundLife", dir = "soundlife"), list(name = "SPAC", dir = "spac")
)
key_types <- c("CD8+ RTE", "CD4+ RTE", "CD8+ T Naive", "CD4+ T Naive",
               "gd T Vd2 GZMK+", "CD8+ Tem GZMK+", "CD8+ Tem GZMB+", "CD8+ Temra")

cmv_coef <- function(co) {  # first CMV row of an lm coefficient matrix, NA if none
  r <- grep("^CMV", rownames(co), value = TRUE)[1]
  if (!is.na(r) && nzchar(r)) r else NA_character_
}

fit_term <- function(fit, term) {
  co <- summary(fit)$coefficients
  if (!term %in% rownames(co)) return(c(NA_real_, NA_real_))
  c(co[term, 1], co[term, ncol(co)])  # beta, p (last col is Pr(>|t|))
}

traj_rows <- list(); tem_rows <- list()
for (co in cohorts) {
  counts_path <- file.path(comp_dir, paste0(co$dir, "_ct_counts.csv"))
  meta_path <- file.path(pb_base, co$dir, "sample_metadata.csv")
  if (!file.exists(counts_path)) { cat("missing:", counts_path, "\n"); next }
  if (!file.exists(meta_path)) { cat("missing:", meta_path, "\n"); next }
  cnt <- read.csv(counts_path, stringsAsFactors = FALSE)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  d <- merge(cnt, meta[, c("donor_id", "age", "sex", "CMV")], by = "donor_id")
  d$age <- as.numeric(d$age); d$age10 <- d$age / 10
  tot <- tapply(d$n_cells, d$donor_id, sum)
  d$total <- tot[d$donor_id]
  d$prop <- d$n_cells / d$total
  d$sex <- factor(d$sex)
  d$CMV <- factor(d$CMV)
  has_cmv <- length(unique(d$CMV)) > 1 && !all(unique(d$CMV) %in% c("No_info", NA))

  for (ct in key_types) {
    s <- d[d$cell_type == ct, ]
    if (nrow(s) < 10) next
    fml <- if (has_cmv) prop ~ age10 + sex + CMV else prop ~ age10 + sex
    fit <- tryCatch(lm(fml, s), error = function(e) NULL)
    if (is.null(fit)) next
    ab <- fit_term(fit, "age10")
    cb <- if (has_cmv) {
      r <- cmv_coef(summary(fit)$coefficients)
      if (is.na(r)) c(NA_real_, NA_real_) else fit_term(fit, r)
    } else c(NA_real_, NA_real_)
    traj_rows[[length(traj_rows) + 1]] <- data.frame(
      cohort = co$name, cell_type = ct, n_donors = nrow(s),
      mean_prop = round(mean(s$prop), 6),
      age_beta_per10y_prop = round(ab[1], 6),
      age_rel_change_per10y = round(ab[1] / mean(s$prop), 3),
      age_p = ab[2], cmv_beta = cb[1], cmv_p = cb[2])
  }

  # Tem GZMK+/GZMB+ ratio
  gzk <- d[d$cell_type == "CD8+ Tem GZMK+", c("donor_id", "n_cells")]
  gzb <- d[d$cell_type == "CD8+ Tem GZMB+", c("donor_id", "n_cells")]
  both <- merge(gzk, gzb, by = "donor_id", suffixes = c("_gzk", "_gzb"))
  both <- merge(both, meta[, c("donor_id", "age", "sex", "CMV")], by = "donor_id")
  if (nrow(both) >= 10) {
    both$age <- as.numeric(both$age); both$age10 <- both$age / 10
    both$sex <- factor(both$sex); both$CMV <- factor(both$CMV)
    both$ratio <- log((both$n_cells_gzk + 0.5) / (both$n_cells_gzb + 0.5))
    has_cmv2 <- length(unique(both$CMV)) > 1 && !all(unique(both$CMV) %in% c("No_info", NA))
    fml2 <- if (has_cmv2) ratio ~ age10 + sex + CMV else ratio ~ age10 + sex
    fit2 <- tryCatch(lm(fml2, both), error = function(e) NULL)
    if (!is.null(fit2)) {
      ab <- fit_term(fit2, "age10")
      cb <- if (has_cmv2) {
        r <- cmv_coef(summary(fit2)$coefficients)
        if (is.na(r)) c(NA_real_, NA_real_) else fit_term(fit2, r)
      } else c(NA_real_, NA_real_)
      tem_rows[[length(tem_rows) + 1]] <- data.frame(
        cohort = co$name, n_donors = nrow(both),
        median_ratio = round(median(exp(both$ratio)), 3),
        age_beta_log_per10y = round(ab[1], 4), age_p = ab[2],
        cmv_beta_log = cb[1], cmv_p = cb[2])
    }
  }
}

traj <- do.call(rbind, traj_rows)
tem <- do.call(rbind, tem_rows)
write.csv(traj, file.path(out_dir, "composition_age_trajectories.csv"), row.names = FALSE)
write.csv(tem, file.path(out_dir, "tem_ratio_age_cmv.csv"), row.names = FALSE)
cat("written:", file.path(out_dir, "composition_age_trajectories.csv"), "\n")
cat("written:", file.path(out_dir, "tem_ratio_age_cmv.csv"), "\n")
print(traj); print(tem)
