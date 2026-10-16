# Association-profile concordance between discovery (ABF300) and query cohorts.
#
# For each cell type, takes the optimal program matching from the projection
# concordance file (proj program <-> refit program), then compares covariate
# association profiles across matched program pairs (wide assoc schema:
# age_beta_per10y_sd, sex_beta_sd, cmv_beta_sd):
#   - Spearman cor of age betas (per cell type and pooled)
#   - sex beta sign agreement
#   - CMV beta concordance (cohorts where CMV is informative)
#
# Usage: Rscript assoc_concordance.R <disc_assoc.csv> <query_assoc.csv> <concordance.csv> <out.csv>

args <- commandArgs(TRUE)
disc_csv <- args[1]
query_csv <- args[2]
conc_csv <- args[3]
out_csv <- args[4]

d <- read.csv(disc_csv, stringsAsFactors = FALSE)
q <- read.csv(query_csv, stringsAsFactors = FALSE)
m <- read.csv(conc_csv, stringsAsFactors = FALSE)

getb <- function(df, ct, prog, col) {
  r <- df[df$cell_type == ct & df$program == prog, ]
  if (!nrow(r) || !col %in% names(r)) return(NA_real_)
  v <- r[[col]][1]
  if (is.na(v) || is.nan(v)) NA_real_ else v
}

has_cmv_d <- any(!is.na(d$cmv_beta_sd))
has_cmv_q <- any(!is.na(q$cmv_beta_sd))

rows <- list()
for (i in seq_len(nrow(m))) {
  ct <- m$cell_type[i]
  p_proj <- m$proj[i]   # discovery program id "<ct>_<k>"
  p_ref <- m$refit[i]   # query refit program id "<ct>_<k>"
  rows[[length(rows) + 1]] <- data.frame(
    cell_type = ct, proj_program = p_proj, refit_program = p_ref,
    concordance_cor = m$cor[i],
    age_disc = getb(d, ct, p_proj, "age_beta_per10y_sd"),
    age_query = getb(q, ct, p_ref, "age_beta_per10y_sd"),
    sex_disc = getb(d, ct, p_proj, "sex_beta_sd"),
    sex_query = getb(q, ct, p_ref, "sex_beta_sd"),
    cmv_disc = getb(d, ct, p_proj, "cmv_beta_sd"),
    cmv_query = getb(q, ct, p_ref, "cmv_beta_sd"))
}
out <- do.call(rbind, rows)

summ <- do.call(rbind, lapply(c("age", "sex", "cmv"), function(tr) {
  dc <- out[[paste0(tr, "_disc")]]; qc <- out[[paste0(tr, "_query")]]
  ok <- !is.na(dc) & !is.na(qc)
  if (sum(ok) < 10) return(data.frame(term = tr, n_pairs = sum(ok), spearman = NA_real_,
                                      pearson = NA_real_, sign_agree = NA_real_))
  data.frame(term = tr, n_pairs = sum(ok),
             spearman = suppressWarnings(cor(dc[ok], qc[ok], method = "spearman")),
             pearson = suppressWarnings(cor(dc[ok], qc[ok])),
             sign_agree = mean(sign(dc[ok]) == sign(qc[ok])))
}))

write.csv(out, out_csv, row.names = FALSE)
cat("assoc concordance written:", out_csv, "\n")
print(summ)
age <- out[!is.na(out$age_disc) & !is.na(out$age_query), ]
if (nrow(age) >= 10) {
  per_ct <- do.call(rbind, lapply(split(age, age$cell_type), function(s) {
    data.frame(cell_type = s$cell_type[1], n = nrow(s),
               spearman = if (nrow(s) >= 5) suppressWarnings(cor(s$age_disc, s$age_query, method = "spearman")) else NA_real_)
  }))
  per_ct <- per_ct[order(-per_ct$spearman), ]
  cat("\nage profile cor by cell type (top 5):\n"); print(head(per_ct, 5))
  cat("\nage profile cor by cell type (bottom 5):\n"); print(tail(per_ct, 5))
}
