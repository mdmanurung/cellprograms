# Discovery-cohort sensitivity of the conservation panel.
#
# Recomputes the tier-1/tier-2 panel metrics with OneK1K (981 donors) as the
# discovery freeze, using the refreeze projection outputs
# (<dir>_vs_onek1k_concordance.csv, <dir>_vs_onek1k_assoc_concordance.csv),
# and joins them with the registered ABF300-discovery panel for comparison.
#
# Usage: Rscript discovery_sensitivity.R <out_csv>

args <- commandArgs(TRUE)
out_csv <- if (length(args)) args[1] else "/workspace/aging_pbmc/discovery_sensitivity.csv"

base <- "/mnt/results/aging_pbmc_repro"
shared <- "/mnt/shared-workspace/aging_pbmc/pseudobulk"
local_pb <- "/workspace/aging_pbmc/pseudobulk"

panel_metrics <- function(pairs) {
  rows <- list()
  for (p in pairs) {
    conc_path <- p$conc
    if (!file.exists(conc_path)) { cat("missing:", conc_path, "\n"); next }
    m <- read.csv(conc_path, stringsAsFactors = FALSE)
    r <- data.frame(query = p$query, query_n_donors = p$query_n,
                    discovery = p$disc, disc_n_donors = p$disc_n,
                    n_matched_programs = nrow(m),
                    median_abs_cor = round(median(abs(m$cor)), 3),
                    frac_abs_cor_ge_0.5 = round(mean(abs(m$cor) >= 0.5), 3),
                    sign_flip_rate = round(mean(m$cor < 0), 3))
    if (file.exists(p$assoc)) {
      o <- read.csv(p$assoc, stringsAsFactors = FALSE)
      dc <- o$age_disc; qc <- o$age_query
      ok <- !is.na(dc) & !is.na(qc)
      if (sum(ok) >= 10) {
        r$assoc_spearman_age <- round(suppressWarnings(cor(dc[ok], qc[ok], method = "spearman")), 3)
        r$assoc_sign_agree_age <- round(mean(sign(dc[ok]) == sign(qc[ok])), 3)
      } else { r$assoc_spearman_age <- NA_real_; r$assoc_sign_agree_age <- NA_real_ }
    }
    rows[[length(rows) + 1]] <- r
  }
  do.call(rbind, rows)
}

# --- ABF300-discovery panel (registered; SoundLife pair is vs KCL) ---
abf_pairs <- list(
  list(query = "SPAC", disc = "ABF300", query_n = 33, disc_n = 166,
       conc = file.path(base, "spac/spac_vs_abf300_concordance.csv"),
       assoc = file.path(base, "spac/spac_assoc_concordance.csv")),
  list(query = "UCSF", disc = "ABF300", query_n = 99, disc_n = 166,
       conc = file.path(base, "ucsf/ucsf_vs_abf300_concordance.csv"),
       assoc = file.path(base, "ucsf/ucsf_assoc_concordance.csv")),
  list(query = "AIDA", disc = "ABF300", query_n = 625, disc_n = 166,
       conc = file.path(base, "aida/aida_vs_abf300_concordance.csv"),
       assoc = file.path(base, "aida/aida_assoc_concordance.csv")),
  list(query = "CIMA", disc = "ABF300", query_n = 421, disc_n = 166,
       conc = file.path(base, "cima/cima_vs_abf300_concordance.csv"),
       assoc = file.path(base, "cima/cima_assoc_concordance.csv")),
  list(query = "OneK1K", disc = "ABF300", query_n = 981, disc_n = 166,
       conc = file.path(shared, "onek1k/onek1k_vs_abf300_concordance.csv"),
       assoc = file.path(shared, "onek1k/onek1k_assoc_concordance.csv")),
  list(query = "SoundLife", disc = "KCL", query_n = 96, disc_n = 188,
       conc = file.path(base, "soundlife/soundlife_vs_kcl_concordance.csv"),
       assoc = file.path(base, "soundlife/soundlife_assoc_concordance.csv"))
)

# --- OneK1K-discovery panel (refreeze; all 7 non-OneK1K cohorts) ---
k1k_pairs <- list(
  list(query = "ABF300", disc = "OneK1K", query_n = 166, disc_n = 981,
       conc = file.path(base, "abf300/abf300_vs_onek1k_concordance.csv"),
       assoc = file.path(base, "abf300/abf300_vs_onek1k_assoc_concordance.csv")),
  list(query = "SPAC", disc = "OneK1K", query_n = 33, disc_n = 981,
       conc = file.path(local_pb, "spac/spac_vs_onek1k_concordance.csv"),
       assoc = file.path(local_pb, "spac/spac_vs_onek1k_assoc_concordance.csv")),
  list(query = "UCSF", disc = "OneK1K", query_n = 99, disc_n = 981,
       conc = file.path(local_pb, "ucsf/ucsf_vs_onek1k_concordance.csv"),
       assoc = file.path(local_pb, "ucsf/ucsf_vs_onek1k_assoc_concordance.csv")),
  list(query = "AIDA", disc = "OneK1K", query_n = 625, disc_n = 981,
       conc = file.path(local_pb, "aida/aida_vs_onek1k_concordance.csv"),
       assoc = file.path(local_pb, "aida/aida_vs_onek1k_assoc_concordance.csv")),
  list(query = "CIMA", disc = "OneK1K", query_n = 421, disc_n = 981,
       conc = file.path(local_pb, "cima/cima_vs_onek1k_concordance.csv"),
       assoc = file.path(local_pb, "cima/cima_vs_onek1k_assoc_concordance.csv")),
  list(query = "KCL", disc = "OneK1K", query_n = 188, disc_n = 981,
       conc = file.path(local_pb, "kcl/kcl_vs_onek1k_concordance.csv"),
       assoc = file.path(local_pb, "kcl/kcl_vs_onek1k_assoc_concordance.csv")),
  list(query = "SoundLife", disc = "OneK1K", query_n = 96, disc_n = 981,
       conc = file.path(local_pb, "soundlife/soundlife_vs_onek1k_concordance.csv"),
       assoc = file.path(local_pb, "soundlife/soundlife_vs_onek1k_assoc_concordance.csv"))
)

panel_abf <- panel_metrics(abf_pairs)
panel_k1k <- panel_metrics(k1k_pairs)

# join per query cohort on shared metrics (KCL/SoundLife have no ABF300-disc
# pair except SoundLife-vs-KCL, kept as its own design)
mets <- c("n_matched_programs", "median_abs_cor", "frac_abs_cor_ge_0.5",
          "sign_flip_rate", "assoc_spearman_age", "assoc_sign_agree_age")
a <- panel_abf[, c("query", mets)]; names(a)[2:length(a)] <- paste0("abf300_", mets)
k <- panel_k1k[, c("query", mets)]; names(k)[2:length(k)] <- paste0("onek1k_", mets)
cmp <- merge(a, k, by = "query", all = TRUE)
cmp <- cmp[order(cmp$query), ]
write.csv(cmp, out_csv, row.names = FALSE)
cat("discovery sensitivity table written:", out_csv, "\n")
print(cmp)
