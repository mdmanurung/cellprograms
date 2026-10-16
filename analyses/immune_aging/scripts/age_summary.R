# Age-association summary across cohorts (refit-side program associations).
#
# Per cohort: programs tested, fraction with age FDR < 0.05, median beta of
# significant hits. Per cell type: cross-cohort sign consistency of the
# strongest (lowest-FDR) program's age beta.
#
# Usage: Rscript age_summary.R <out_csv>

cohorts <- list(
  list(name = "ABF300", n = 166, csv = "/workspace/aging_pbmc/pseudobulk/abf300/abf300_full_assoc.csv"),
  list(name = "OneK1K", n = 981, csv = "/mnt/shared-workspace/aging_pbmc/pseudobulk/onek1k/onek1k_full_assoc.csv"),
  list(name = "AIDA", n = 625, csv = "/workspace/aging_pbmc/pseudobulk/aida/aida_full_assoc.csv"),
  list(name = "CIMA", n = 421, csv = "/workspace/aging_pbmc/pseudobulk/cima/cima_full_assoc.csv"),
  list(name = "KCL", n = 188, csv = "/workspace/aging_pbmc/pseudobulk/kcl/kcl_full_assoc.csv"),
  list(name = "UCSF", n = 99, csv = "/workspace/aging_pbmc/pseudobulk/ucsf/ucsf_full_assoc.csv"),
  list(name = "SoundLife", n = 96, csv = "/workspace/aging_pbmc/pseudobulk/soundlife/soundlife_full_assoc.csv"),
  list(name = "SPAC", n = 33, csv = "/workspace/aging_pbmc/pseudobulk/spac/spac_full_assoc.csv")
)

out <- list()
per_ct_rows <- list()
for (co in cohorts) {
  if (!file.exists(co$csv)) { cat("missing:", co$csv, "\n"); next }
  d <- read.csv(co$csv, stringsAsFactors = FALSE)
  sig <- d[!is.na(d$age_p_fdr) & d$age_p_fdr < 0.05, ]
  out[[length(out) + 1]] <- data.frame(
    cohort = co$name, n_donors = co$n, n_programs = nrow(d),
    n_age_fdr05 = nrow(sig), frac_age_fdr05 = round(nrow(sig) / max(nrow(d), 1), 3),
    median_beta_sig = if (nrow(sig)) round(median(sig$age_beta_per10y_sd), 3) else NA_real_,
    n_neg_sig = sum(sig$age_beta_per10y_sd < 0), n_pos_sig = sum(sig$age_beta_per10y_sd > 0))
  # strongest program per cell type
  for (ct in unique(d$cell_type)) {
    s <- d[d$cell_type == ct & !is.na(d$age_p_fdr), ]
    if (!nrow(s)) next
    w <- s[which.min(s$age_p_fdr), ]
    per_ct_rows[[length(per_ct_rows) + 1]] <- data.frame(
      cohort = co$name, cell_type = ct, program = w$program,
      age_beta = round(w$age_beta_per10y_sd, 3), age_p_fdr = w$age_p_fdr)
  }
}
tab <- do.call(rbind, out)
write.csv(tab, "/workspace/aging_pbmc/age_by_cohort.csv", row.names = FALSE)
cat("\n=== Age associations by cohort (refit-side) ===\n")
print(tab, row.names = FALSE)

pc <- do.call(rbind, per_ct_rows)
# cross-cohort sign consistency per cell type (cohorts with FDR<0.05 only)
sig <- pc[pc$age_p_fdr < 0.05, ]
cons <- do.call(rbind, lapply(split(sig, sig$cell_type), function(s) {
  data.frame(cell_type = s$cell_type[1], n_cohorts_sig = nrow(s),
             n_cohorts_total = length(unique(pc$cell_type[pc$cell_type == s$cell_type[1]])),
             n_neg = sum(s$age_beta < 0), n_pos = sum(s$age_beta > 0),
             median_beta = round(median(s$age_beta), 3),
             cohorts = paste(s$cohort[order(s$age_beta)], collapse = ","))
}))
cons <- cons[order(-cons$n_cohorts_sig, cons$cell_type), ]
write.csv(cons, "/workspace/aging_pbmc/age_by_celltype_consistency.csv", row.names = FALSE)
cat("\n=== Cell types with age-associated programs in >=2 cohorts (FDR<0.05) ===\n")
print(cons[cons$n_cohorts_sig >= 2, ], row.names = FALSE)
