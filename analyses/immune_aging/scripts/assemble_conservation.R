# Assemble the cross-cohort conservation table (3-tier panel).
#
# Tier 1 (co-primary): MR-1 score-space concordance (matched |cor| between
#   projected and refit program scores) + association-profile replication
#   (Spearman cor of covariate betas across matched programs).
# Tier 2: descriptive matched-concordance distribution stats (sign flips,
#   fraction >= 0.5) — principal-angle/subspace metrics deferred.
# Tier 3: gene-set overlap — deferred to the report stage (needs frozen
#   top-gene sets per program; computed separately).
#
# Usage: Rscript assemble_conservation.R <out_dir>
# Expects per-pair outputs under /workspace/aging_pbmc/pseudobulk/<query>/:
#   <query>_vs_<disc>_concordance.csv and <query>_assoc_concordance.csv

args <- commandArgs(TRUE)
out_dir <- if (length(args)) args[1] else "/workspace/aging_pbmc/conservation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

base <- "/workspace/aging_pbmc/pseudobulk"
pairs <- list(
  data.frame(query = "SPAC", disc = "ABF300", query_n = 33, disc_n = 166,
             conc = "spac/spac_vs_abf300_concordance.csv", assoc = "spac/spac_assoc_concordance.csv"),
  data.frame(query = "UCSF", disc = "ABF300", query_n = 99, disc_n = 166,
             conc = "ucsf/ucsf_vs_abf300_concordance.csv", assoc = "ucsf/ucsf_assoc_concordance.csv"),
  data.frame(query = "AIDA", disc = "ABF300", query_n = 625, disc_n = 166,
             conc = "aida/aida_vs_abf300_concordance.csv", assoc = "aida/aida_assoc_concordance.csv"),
  data.frame(query = "CIMA", disc = "ABF300", query_n = 421, disc_n = 166,
             conc = "cima/cima_vs_abf300_concordance.csv", assoc = "cima/cima_assoc_concordance.csv"),
  data.frame(query = "OneK1K", disc = "ABF300", query_n = 981, disc_n = 166,
             conc = "/mnt/shared-workspace/aging_pbmc/pseudobulk/onek1k/onek1k_vs_abf300_concordance.csv",
             assoc = "/mnt/shared-workspace/aging_pbmc/pseudobulk/onek1k/onek1k_assoc_concordance.csv"),
  data.frame(query = "SoundLife", disc = "KCL", query_n = 96, disc_n = 188,
             conc = "soundlife/soundlife_vs_kcl_concordance.csv", assoc = "soundlife/soundlife_assoc_concordance.csv")
)

rows <- list()
for (p in pairs) {
  conc_path <- file.path(base, p$conc)
  assoc_path <- file.path(base, p$assoc)
  if (!file.exists(conc_path)) next
  m <- read.csv(conc_path, stringsAsFactors = FALSE)
  r <- data.frame(query = p$query, query_n_donors = p$query_n,
                  discovery = p$disc, disc_n_donors = p$disc_n,
                  n_matched_programs = nrow(m),
                  median_abs_cor = round(median(abs(m$cor)), 3),
                  frac_abs_cor_ge_0.5 = round(mean(abs(m$cor) >= 0.5), 3),
                  sign_flip_rate = round(mean(m$cor < 0), 3))
  if (file.exists(assoc_path)) {
    o <- read.csv(assoc_path, stringsAsFactors = FALSE)
    for (tm in c("age", "sex", "cmv")) {
      dc <- o[[paste0(tm, "_disc")]]; qc <- o[[paste0(tm, "_query")]]
      ok <- !is.na(dc) & !is.na(qc)
      if (sum(ok) >= 10) {
        r[[paste0("assoc_spearman_", tm)]] <- round(suppressWarnings(cor(dc[ok], qc[ok], method = "spearman")), 3)
        r[[paste0("assoc_sign_agree_", tm)]] <- round(mean(sign(dc[ok]) == sign(qc[ok])), 3)
      } else {
        r[[paste0("assoc_spearman_", tm)]] <- NA_real_
        r[[paste0("assoc_sign_agree_", tm)]] <- NA_real_
      }
    }
  }
  rows[[length(rows) + 1]] <- r
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "conservation_panel.csv"), row.names = FALSE)
cat("conservation panel written:", file.path(out_dir, "conservation_panel.csv"), "\n")
print(tab)
