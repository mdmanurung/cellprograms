#' cellprograms: Cell-Type-Specific Transcriptional Programs via EBMF
#'
#' Discover cell-type-specific transcriptional programs from pseudobulk
#' single-cell RNA-seq data using empirical-Bayes matrix factorization
#' ('flashier'), with canonicalization, stability assessment, cross-cell-type
#' coordination detection, shared-vs-specific program assessment ('mashr'),
#' and artifact-gene masking. Programs are discovered locally within each cell
#' type; cross-cell-type coordination is investigated subsequently, never
#' assumed.
#'
#' @keywords internal
#' @import stats utils
"_PACKAGE"

## ggplot2 data-column references in R/plots.R (NSE; not global variables)
utils::globalVariables(c(
  "a", "ambient", "b", "cell_type", "cmp_lab", "direction", "donor",
  "estimate", "fraction", "gene", "lab", "label", "loading",
  "median_score_cor", "min_angle_deg", "n_shared_05", "node_stat", "p",
  "pair", "program", "pve", "recovery_freq", "score", "se", "set",
  "sign_consistency", "sig", "stat", "to", "top_gene_jaccard",
  "underpowered", "value", "weight", "x", "x0", "x1", "xmax", "xmin",
  "y", "y0", "y1", "ymax", "ymax2", "ymin", "ymin2", "level", "n",
  "programs"
))

