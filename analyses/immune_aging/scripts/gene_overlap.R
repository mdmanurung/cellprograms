# Tier-3 gene-set overlap: top-N gene Jaccard for matched program pairs.
#
# For each pair (discovery fit, query refit fit) and each matched program pair
# from the projection concordance file: Jaccard of top-20 genes by |loading|.
#
# Usage: Rscript gene_overlap.R <out_csv>

suppressPackageStartupMessages(library(cellprograms))

base <- "/workspace/aging_pbmc/pseudobulk"
pairs <- list(
  list(query = "SPAC", qdir = "spac", disc = "abf300", conc = "spac/spac_vs_abf300_concordance.csv"),
  list(query = "UCSF", qdir = "ucsf", disc = "abf300", conc = "ucsf/ucsf_vs_abf300_concordance.csv"),
  list(query = "AIDA", qdir = "aida", disc = "abf300", conc = "aida/aida_vs_abf300_concordance.csv"),
  list(query = "CIMA", qdir = "cima", disc = "abf300", conc = "cima/cima_vs_abf300_concordance.csv"),
  list(query = "SoundLife", qdir = "soundlife", disc = "kcl", conc = "soundlife/soundlife_vs_kcl_concordance.csv")
)
N <- 20

rows <- list()
for (p in pairs) {
  conc_path <- file.path(base, p$conc)
  if (!file.exists(conc_path)) next
  m <- read.csv(conc_path, stringsAsFactors = FALSE)
  disc <- readRDS(file.path(base, p$disc, paste0(p$disc, "_full_fit.rds")))
  qfit <- readRDS(file.path(base, p$qdir, paste0(tolower(p$query), "_full_fit.rds")))
  for (i in seq_len(nrow(m))) {
    ct <- m$cell_type[i]
    if (!ct %in% names(disc$fits) || !ct %in% names(qfit$fits)) next
    Wd <- disc$fits[[ct]]$flash$F_pm
    Wq <- qfit$fits[[ct]]$flash$F_pm
    kd <- suppressWarnings(as.integer(sub(".*_", "", m$proj[i])))
    kq <- suppressWarnings(as.integer(sub(".*_", "", m$refit[i])))
    if (is.na(kd) || is.na(kq) || kd > ncol(Wd) || kq > ncol(Wq)) next
    td <- head(rownames(Wd)[order(-abs(Wd[, kd]))], N)
    tq <- head(rownames(Wq)[order(-abs(Wq[, kq]))], N)
    rows[[length(rows) + 1]] <- data.frame(
      query = p$query, cell_type = ct, proj = m$proj[i], refit = m$refit[i],
      score_cor = m$cor[i], jaccard_top20 = length(intersect(td, tq)) / length(union(td, tq)))
  }
}
out <- do.call(rbind, rows)
write.csv(out, "/workspace/aging_pbmc/gene_overlap_top20.csv", row.names = FALSE)
cat("gene overlap written: /workspace/aging_pbmc/gene_overlap_top20.csv\n")
summ <- do.call(rbind, lapply(split(out, out$query), function(s) {
  data.frame(query = s$query[1], n = nrow(s),
             median_jaccard = round(median(s$jaccard_top20), 3),
             frac_jaccard_ge_0.2 = round(mean(s$jaccard_top20 >= 0.2), 3))
}))
print(summ, row.names = FALSE)
