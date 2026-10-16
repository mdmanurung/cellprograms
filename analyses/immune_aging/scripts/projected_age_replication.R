# Projected-score age replication: test discovery programs' age associations
# directly on their PROJECTED scores in each query cohort.
#
# For each query cohort with a saved _vs_abf300_activity.rds:
#   - per cell type and discovery program: lm(projected_score ~ age + sex)
#   - compare query-side age beta with the discovery-side (ABF300 refit) beta
#   - report Spearman cor across programs per ct and pooled, sign agreement,
#     and replication of significant discovery hits (FDR<0.05) in the query
#
# Usage: Rscript projected_age_replication.R <out_csv>

suppressPackageStartupMessages(library(cellprograms))

base <- "/workspace/aging_pbmc/pseudobulk"
disc_assoc <- read.csv(file.path(base, "abf300/abf300_full_assoc.csv"), stringsAsFactors = FALSE)
queries <- list(
  list(name = "SPAC", n = 33, dir = "spac", tag = "spac_vs_abf300"),
  list(name = "UCSF", n = 99, dir = "ucsf", tag = "ucsf_vs_abf300"),
  list(name = "AIDA", n = 625, dir = "aida", tag = "aida_vs_abf300"),
  list(name = "CIMA", n = 421, dir = "cima", tag = "cima_vs_abf300"),
  list(name = "OneK1K", n = 981, dir = "onek1k", tag = "onek1k_vs_abf300")
)
# OneK1K projection outputs live on shared storage
act_path_fix <- function(q) {
  p <- file.path(base, q$dir, paste0(q$tag, "_activity.rds"))
  if (file.exists(p)) return(p)
  p2 <- file.path("/mnt/shared-workspace/aging_pbmc/pseudobulk", q$dir, paste0(q$tag, "_activity.rds"))
  if (file.exists(p2)) return(p2)
  p
}

all_rows <- list()
summ <- list()
for (q in queries) {
  act_path <- act_path_fix(q)
  if (!file.exists(act_path)) { cat("missing:", act_path, "\n"); next }
  act <- readRDS(act_path)
  ids <- rownames(act)
  dvec <- sub("\\|[^|]*$", "", ids)
  cts_vec <- sub("^[^|]*\\|", "", ids)
  meta <- read.csv(file.path(base, q$dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
  age_map <- stats::setNames(meta$age, meta$donor_id)
  sex_map <- stats::setNames(meta$sex, meta$donor_id)
  for (ct in unique(cts_vec)) {
    ids_ct <- ids[cts_vec == ct]
    sub <- act[ids_ct, , drop = FALSE]
    # drop all-NA columns (programs not active in this ct)
    keep <- colSums(!is.na(sub)) > nrow(sub) / 2
    sub <- sub[, keep, drop = FALSE]
    if (!ncol(sub)) next
    ag <- as.numeric(age_map[sub("\\|[^|]*$", "", ids_ct)])
    sx <- sex_map[sub("\\|[^|]*$", "", ids_ct)]
    ok <- !is.na(ag)
    if (sum(ok) < 20) next
    sx <- factor(sx)
    if (length(levels(sx)) < 2) sx <- NULL
    for (p in colnames(sub)) {
      y <- sub[ok, p]
      if (sum(!is.na(y)) < 20) next
      dd <- if (is.null(sx)) data.frame(score = y, age = ag[ok]) else data.frame(score = y, age = ag[ok], sex = sx)
      fml <- if (is.null(sx)) score ~ age else score ~ age + sex
      fit_lm <- tryCatch(lm(fml, dd, na.action = na.omit), error = function(e) NULL)
      if (is.null(fit_lm)) next
      co <- summary(fit_lm)$coefficients
      if (!"age" %in% rownames(co)) next
      all_rows[[length(all_rows) + 1]] <- data.frame(
        query = q$name, cell_type = ct, program = p,
        age_beta_query_proj = co["age", 1], age_p_query = co["age", 4])
    }
  }
}
out <- do.call(rbind, all_rows)
# merge with discovery-side betas
out$age_beta_disc <- disc_assoc$age_beta_per10y_sd[match(out$program, disc_assoc$program)]
out <- out[!is.na(out$age_beta_disc), ]
write.csv(out, "/workspace/aging_pbmc/projected_age_replication.csv", row.names = FALSE)

# summary per query
for (qn in unique(out$query)) {
  s <- out[out$query == qn, ]
  sp <- suppressWarnings(cor(s$age_beta_disc, s$age_beta_query_proj, method = "spearman"))
  # replication of significant discovery hits
  disc_sig <- s[abs(s$age_beta_disc) > 0 & !is.na(s$age_p_query), ]
  # discovery significance from assoc file
  dsig <- disc_assoc$age_p_fdr[match(s$program, disc_assoc$program)]
  hits <- s[!is.na(dsig) & dsig < 0.05, ]
  repl <- if (nrow(hits)) mean(hits$age_p_query < 0.05 & sign(hits$age_beta_query_proj) == sign(hits$age_beta_disc)) else NA
  cat(qn, ": n =", nrow(s), "| spearman =", round(sp, 3),
      "| disc-sig hits:", nrow(hits), "| replicated in query (p<0.05, same sign):", round(repl, 3), "\n")
}
