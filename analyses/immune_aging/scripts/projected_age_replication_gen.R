# Projected-score age replication, parameterized by discovery freeze.
#
# For each query cohort with a saved <dir>_vs_<disc>_activity.rds:
#   - per cell type and discovery program: lm(projected_score ~ age + sex)
#   - compare query-side age beta with the discovery-side beta
#   - report Spearman cor across programs, and replication of significant
#     discovery hits (FDR<0.05) in the query (p<0.05, same sign)
#
# Usage: Rscript projected_age_replication_gen.R <disc> <out_csv>
#   <disc> in {abf300, onek1k}

suppressPackageStartupMessages(library(cellprograms))

args <- commandArgs(trailingOnly = TRUE)
disc <- if (length(args) >= 1) tolower(args[1]) else "abf300"
out_csv <- if (length(args) >= 2) args[2] else
  sprintf("/workspace/aging_pbmc/projected_age_replication_%s_disc.csv", disc)

base <- "/workspace/aging_pbmc/pseudobulk"
shared <- "/mnt/shared-workspace/aging_pbmc/pseudobulk"

# discovery-side association file (wide schema)
disc_assoc_path <- file.path(base, disc, paste0(disc, "_full_assoc.csv"))
if (!file.exists(disc_assoc_path)) {
  disc_assoc_path <- file.path(shared, disc, paste0(disc, "_full_assoc.csv"))
}
stopifnot(file.exists(disc_assoc_path))
disc_assoc <- read.csv(disc_assoc_path, stringsAsFactors = FALSE)

# all cohorts except the discovery itself
all_cohorts <- list(
  list(name = "ABF300", n = 166, dir = "abf300"),
  list(name = "OneK1K", n = 981, dir = "onek1k"),
  list(name = "AIDA", n = 625, dir = "aida"),
  list(name = "CIMA", n = 421, dir = "cima"),
  list(name = "KCL", n = 188, dir = "kcl"),
  list(name = "UCSF", n = 99, dir = "ucsf"),
  list(name = "SoundLife", n = 96, dir = "soundlife"),
  list(name = "SPAC", n = 33, dir = "spac")
)
queries <- Filter(function(q) tolower(q$name) != disc, all_cohorts)

act_path_fix <- function(q) {
  tag <- sprintf("%s_vs_%s", q$dir, disc)
  p <- file.path(base, q$dir, paste0(tag, "_activity.rds"))
  if (file.exists(p)) return(p)
  p2 <- file.path(shared, q$dir, paste0(tag, "_activity.rds"))
  if (file.exists(p2)) return(p2)
  p
}

all_rows <- list()
for (q in queries) {
  act_path <- act_path_fix(q)
  if (!file.exists(act_path)) { cat("missing:", act_path, "\n"); next }
  act <- readRDS(act_path)
  ids <- rownames(act)
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
        discovery = disc, query = q$name, cell_type = ct, program = p,
        age_beta_query_proj = co["age", 1], age_p_query = co["age", 4])
    }
  }
}
out <- do.call(rbind, all_rows)
# merge with discovery-side betas
out$age_beta_disc <- disc_assoc$age_beta_per10y_sd[match(out$program, disc_assoc$program)]
out <- out[!is.na(out$age_beta_disc), ]
write.csv(out, out_csv, row.names = FALSE)

# summary per query
cat("=== projected age replication, discovery =", disc, "===\n")
for (qn in unique(out$query)) {
  s <- out[out$query == qn, ]
  sp <- suppressWarnings(cor(s$age_beta_disc, s$age_beta_query_proj, method = "spearman"))
  dsig <- disc_assoc$age_p_fdr[match(s$program, disc_assoc$program)]
  hits <- s[!is.na(dsig) & dsig < 0.05, ]
  repl <- if (nrow(hits)) mean(hits$age_p_query < 0.05 & sign(hits$age_beta_query_proj) == sign(hits$age_beta_disc)) else NA
  cat(qn, ": n =", nrow(s), "| spearman =", round(sp, 3),
      "| disc-sig hits:", nrow(hits), "| replicated in query (p<0.05, same sign):", round(repl, 3), "\n")
}
cat("written:", out_csv, "\n")
