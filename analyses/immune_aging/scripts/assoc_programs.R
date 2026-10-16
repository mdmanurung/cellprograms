# Association analysis: program scores vs age, sex, CMV per cell type.
#
# For each cell type and program: lm(score ~ age + sex + CMV) over donors with
# complete covariates; age beta reported per 10 years in SD units of the score,
# sex and CMV betas in SD units; BH-FDR within cell type across programs.
#
# Usage: Rscript assoc_programs.R <fit.rds> <pseudobulk_dir> <out_csv>
suppressPackageStartupMessages(library(cellprograms))

args <- commandArgs(TRUE)
fit_path <- args[1]
pb_dir <- args[2]
out_csv <- args[3]

fit <- readRDS(fit_path)
meta <- read.csv(file.path(pb_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
if (!"observation_id" %in% names(meta)) meta$observation_id <- meta$donor_id

rows <- list()
for (ct in names(fit$scores)) {
  S <- fit$scores[[ct]]
  if (is.null(S) || nrow(S) < 10) next
  d <- data.frame(observation_id = rownames(S), stringsAsFactors = FALSE)
  m <- meta[match(d$observation_id, meta$observation_id), , drop = FALSE]
  d$age <- m$age
  d$sex <- m$sex
  d$CMV <- m$CMV
  d <- d[complete.cases(d$age, d$sex), ]
  if (nrow(d) < 10) next
  d$sex <- factor(d$sex)
  has_cmv <- !all(d$CMV %in% c("No_info", NA)) && length(unique(d$CMV[!is.na(d$CMV)])) > 1
  if (has_cmv) d$CMV <- factor(d$CMV)

  for (p in colnames(S)) {
    dd <- d
    dd$score <- S[dd$observation_id, p]
    fit_lm <- tryCatch(
      if (has_cmv) lm(score ~ age + sex + CMV, dd) else lm(score ~ age + sex, dd),
      error = function(e) NULL)
    if (is.null(fit_lm)) next  # degenerate design (e.g. single sex after complete cases)
    co <- summary(fit_lm)$coefficients
    sd_y <- sd(dd$score)
    if (is.na(sd_y) || sd_y == 0) next
    age_beta <- if ("age" %in% rownames(co)) co["age", 1] * 10 / sd_y else NA
    age_p <- if ("age" %in% rownames(co)) co["age", 4] else NA
    sex_beta <- sex_p <- NA
    sl <- grep("^sex", rownames(co))
    if (length(sl)) { sex_beta <- co[sl[1], 1] / sd_y; sex_p <- co[sl[1], 4] }
    cmv_beta <- cmv_p <- NA
    cl <- grep("^CMV", rownames(co))
    if (length(cl)) { cmv_beta <- co[cl[1], 1] / sd_y; cmv_p <- co[cl[1], 4] }
    rows[[length(rows) + 1]] <- data.frame(
      cell_type = ct, program = p, n_donors = nrow(dd),
      age_beta_per10y_sd = age_beta, age_p = age_p,
      sex_beta_sd = sex_beta, sex_p = sex_p,
      cmv_beta_sd = cmv_beta, cmv_p = cmv_p)
  }
}

res <- do.call(rbind, rows)
# BH-FDR within cell type (age association is the primary screen)
res$age_p_fdr <- NA_real_
for (ct in unique(res$cell_type)) {
  i <- res$cell_type == ct & !is.na(res$age_p)
  res$age_p_fdr[i] <- p.adjust(res$age_p[i], method = "BH")
}
write.csv(res, out_csv, row.names = FALSE)
cat("associations written:", out_csv, nrow(res), "rows\n")
