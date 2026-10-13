#!/usr/bin/env Rscript
# Pairing sensitivity for the DIVAS case study (plan step 4b).
#
# The cohort has 2 timepoints (T1/T2) per patient (57 patients, 114 samples).
# The marginal coordination treats samples as independent, so we run two
# checks:
#   (1) T1-only coordination: rerun find_coordinated_programs on the 57
#       baseline samples only -- removes repeated-measures structure entirely.
#   (2) Patient-blocked permutation null: for each observed edge, test its
#       Spearman correlation against a null that permutes whole patients
#       (both timepoints together), preserving within-patient correlation.
#
# Usage: Rscript 07_sensitivity_pairing.R <data_dir> <out_dir>

suppressPackageStartupMessages({
  if (requireNamespace("cellprograms", quietly = TRUE)) {
    library(cellprograms)
  } else {
    pkgload::load_all(Sys.getenv("CELLPROGRAMS_ROOT", unset = "../../.."), quiet = TRUE)
  }
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir <- if (length(args) >= 2) args[2] else "results"
set.seed(1)

fit <- readRDS(file.path(out_dir, "divas6_fit.rds"))
cmp_full <- readRDS(file.path(out_dir, "cmp_marginal.rds"))
meta <- read.csv(file.path(data_dir, "input", "final_sample_info.csv"), stringsAsFactors = FALSE)

# ================= 1. T1-only coordination =================
cat("== T1-only coordination (57 patients, n_perm=999)\n")
t1_ids <- meta$sample_id[meta$timepoint == "T1"]
stopifnot(all(t1_ids %in% fit$observation_ids))
keep <- match(t1_ids, fit$observation_ids)

fit_t1 <- fit
fit_t1$observation_ids <- fit$observation_ids[keep]
for (ct in fit$cell_types) fit_t1$scores[[ct]] <- fit$scores[[ct]][keep, , drop = FALSE]

cmp_t1 <- find_coordinated_programs(fit_t1, n_perm = 999, seed = 1, min_donors = 20)
saveRDS(cmp_t1, file.path(out_dir, "cmp_t1only.rds"))

# -- compare headline calls: full-cohort edges reproduced at T1? --
edge_key <- function(e) paste(e$ct_a, e$prog_a, e$ct_b, e$prog_b, sep = "|")
full_edges <- cmp_full$edges
t1_edges <- cmp_t1$edges
key_full <- vapply(seq_len(nrow(full_edges)), function(i)
  edge_key(full_edges[i, ]), "")
key_t1 <- vapply(seq_len(nrow(t1_edges)), function(i)
  edge_key(t1_edges[i, ]), "")
reproduced <- key_full %in% key_t1
edge_cmp <- data.frame(
  ct_a = full_edges$ct_a, prog_a = full_edges$prog_a,
  ct_b = full_edges$ct_b, prog_b = full_edges$prog_b,
  reproduced_at_T1 = reproduced, stringsAsFactors = FALSE)
write.csv(edge_cmp, file.path(out_dir, "sensitivity_edges_t1.csv"), row.names = FALSE)
cmp_tab_full <- table(cmp_full$edges$cmp)
cmp_tab_t1 <- table(cmp_t1$edges$cmp)
cat(sprintf("edges: full=%d | T1-only=%d | full edges reproduced at T1: %d/%d (%.0f%%)\n",
            nrow(full_edges), nrow(t1_edges), sum(reproduced), length(reproduced),
            100 * mean(reproduced)))
cat(sprintf("CMPs: full=%d | T1-only=%d\n", length(cmp_full$cmps), length(cmp_t1$cmps)))

# ================= 2. Patient-blocked permutation null =================
cat("== patient-blocked permutation null on edge correlations\n")
# observed Spearman correlation for each full-cohort edge
score_of <- function(ct, prog) fit$scores[[ct]][, prog]
pid <- meta$clinical_id[match(fit$observation_ids, meta$sample_id)]
stopifnot(!anyNA(pid))
# patient units: all samples of a patient move together under permutation
# (29 patients contribute a T1/T2 pair, 56 a single timepoint)
pid <- factor(pid)
plevels <- levels(pid)
obs_cor <- vapply(seq_len(nrow(full_edges)), function(i) {
  suppressWarnings(cor(score_of(full_edges$ct_a[i], full_edges$prog_a[i]),
                       score_of(full_edges$ct_b[i], full_edges$prog_b[i]),
                       method = "spearman"))
}, 0)

n_perm <- 999
null_cor <- matrix(NA_real_, nrow = n_perm, ncol = nrow(full_edges))
# precompute per-block patient-level score means (T1+T2 rows move together)
block_scores <- lapply(full_edges$ct_a, function(ct) fit$scores[[ct]])  # per-edge block a
# permutation: shuffle patient labels of block b's rows, in patient units
perm_pvals <- numeric(nrow(full_edges))
set.seed(1)
# build patient-index list once (keyed by clinical_id, matching plevels)
rows_by_patient <- split(seq_along(pid), pid)

for (p in seq_len(n_perm)) {
  # permute patient units for the b-side scores of every edge (one shared
  # permutation per iteration keeps edges comparable)
  new_order <- sample(plevels); names(new_order) <- plevels
  row_perm <- unlist(rows_by_patient[new_order], use.names = FALSE)
  for (j in seq_len(nrow(full_edges))) {
    sb <- fit$scores[[full_edges$ct_b[j]]][, full_edges$prog_b[j]]
    null_cor[p, j] <- suppressWarnings(
      cor(score_of(full_edges$ct_a[j], full_edges$prog_a[j]),
          sb[row_perm], method = "spearman"))
  }
}
for (j in seq_len(nrow(full_edges))) {
  perm_pvals[j] <- (1 + sum(abs(null_cor[, j]) >= abs(obs_cor[j]))) / (n_perm + 1)
}
perm_tab <- data.frame(
  ct_a = full_edges$ct_a, prog_a = full_edges$prog_a,
  ct_b = full_edges$ct_b, prog_b = full_edges$prog_b,
  obs_spearman = obs_cor, p_patient_blocked = perm_pvals,
  stringsAsFactors = FALSE)
write.csv(perm_tab, file.path(out_dir, "sensitivity_patient_blocked_perm.csv"), row.names = FALSE)
cat(sprintf("edges significant at patient-blocked p<0.05: %d/%d (%.0f%%)\n",
            sum(perm_pvals < 0.05), length(perm_pvals), 100 * mean(perm_pvals < 0.05)))
cat("DONE\n")
