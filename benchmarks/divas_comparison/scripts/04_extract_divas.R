# 04_extract_divas.R -- extract DIVAS component scores/loadings/taxonomy
#
# The precomputed results object (Zenodo 17430294) contains:
#   Scores         : samples x 91 component scores (column names encode the
#                    block-membership pattern, e.g. "CD4+CD8+CD14+NK+pro+metab-1")
#   scoresList     : per-block views of the same components
#   Loadings       : per-block feature x component loading matrices
#   jointBasisMapRaw : per-pattern joint basis (samples x k)
#   keyIdxMap / keymapname : pattern bitmask -> component column indices / names
#
# Outputs (CSV, small): component annotation, score matrix, top-100 loadings
# per component per block, pattern table, joint bases.
#
# Usage: Rscript 04_extract_divas.R <data_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir  <- if (length(args) >= 2) args[2] else "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(data_dir, "divas_results_combined_6block_renamed.rds"))
stopifnot(all(c("Scores", "Loadings", "jointBasisMapRaw") %in% names(obj)))

# -- component annotation from column names --
comp_names <- colnames(obj$Scores)
parse_comp <- function(nm) {
  parts <- strsplit(nm, "-", fixed = TRUE)[[1]]
  blocks <- strsplit(parts[1], "+", fixed = TRUE)[[1]]
  data.frame(component = nm,
             pattern = parts[1],
             pattern_index = as.integer(parts[2]),
             n_blocks = length(blocks),
             blocks = paste(blocks, collapse = "+"),
             stringsAsFactors = FALSE)
}
annot <- do.call(rbind, lapply(comp_names, parse_comp))
annot$in_gene_panels <- vapply(strsplit(annot$blocks, "+", fixed = TRUE),
                               function(b) all(c("CD4", "CD8", "CD14", "NK") %in% b) &&
                                 length(b) == sum(c("CD4", "CD8", "CD14", "NK", "pro", "metab") %in% b),
                               logical(1))
# DIVAS-style taxonomy by pattern breadth
annot$divas_class <- ifelse(annot$n_blocks == 6, "joint",
                     ifelse(annot$n_blocks == 1, "individual", "partial"))
cat("components:", nrow(annot), "| patterns:", length(unique(annot$pattern)), "\n")
print(table(annot$divas_class))

write.csv(annot, file.path(out_dir, "divas_components.csv"), row.names = FALSE)

# -- scores (samples x components) --
scores <- obj$Scores
stopifnot(nrow(scores) == 114L, ncol(scores) == nrow(annot))
write.csv(data.frame(sample_id = rownames(scores), scores, check.names = FALSE),
          file.path(out_dir, "divas_scores.csv"), row.names = FALSE)

# -- top-100 |loadings| per component, per block (gene blocks only for the
#    overlap analysis; omics blocks included for completeness) --
top_list <- list()
for (blk in names(obj$Loadings)) {
  L <- obj$Loadings[[blk]]
  if (is.null(L) || ncol(L) == 0L) next
  for (j in seq_len(ncol(L))) {
    o <- order(-abs(L[, j]))[seq_len(min(100L, nrow(L)))]
    top_list[[paste(blk, j, sep = "_")]] <-
      data.frame(block = blk, component = colnames(L)[j], rank = seq_along(o),
                 feature = rownames(L)[o], loading = L[o, j],
                 stringsAsFactors = FALSE)
  }
}
top <- do.call(rbind, top_list)
write.csv(top, file.path(out_dir, "divas_top100_loadings.csv"), row.names = FALSE)

# -- pattern table + joint bases --
pat <- annot[!duplicated(annot$pattern), c("pattern", "n_blocks", "blocks", "divas_class")]
pat$n_components <- as.integer(table(annot$pattern)[pat$pattern])
write.csv(pat, file.path(out_dir, "divas_patterns.csv"), row.names = FALSE)

jb <- obj$jointBasisMapRaw
jb_out <- list()
score_sample_ids <- rownames(obj$Scores)   # joint bases carry no rownames
for (nm in names(jb)) {
  M <- jb[[nm]]
  if (is.null(M) || ncol(M) == 0L) next
  stopifnot(nrow(M) == length(score_sample_ids))
  for (j in seq_len(ncol(M))) {
    jb_out[[paste(nm, j, sep = "_")]] <-
      data.frame(pattern = nm, basis_index = j, sample_id = score_sample_ids,
                 score = M[, j], stringsAsFactors = FALSE)
  }
}
jb_df <- do.call(rbind, jb_out)
write.csv(jb_df, file.path(out_dir, "divas_joint_bases.csv"), row.names = FALSE)

# -- full loadings kept as rds for the overlap analysis --
saveRDS(obj$Loadings, file.path(out_dir, "divas_loadings.rds"))
cat("DONE: extracted", nrow(annot), "components,", nrow(pat), "patterns\n")
