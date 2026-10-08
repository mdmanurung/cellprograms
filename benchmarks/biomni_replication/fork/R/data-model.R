## data-model.R — input data contract
##
## Canonical input: named list of pseudobulk matrices (observations x genes),
## plus optional sample metadata. Rows are biological observations
## (typically participant x timepoint); columns are genes.

#' Create a cell_program_data object
#'
#' @param pseudobulk Named list of numeric matrices (observations x genes).
#'   List names are cell-type labels. Row names must be observation IDs;
#'   column names must be gene IDs.
#' @param sample_metadata Optional data.frame with an `observation_id` column.
#'   Outcome variables may be included but are never used in factorization.
#' @param features Feature-selection mode: `"all"`, `"variable"` (drop
#'   zero-variance genes), or a named list of gene vectors per cell type.
#' @return An object of class `cell_program_data`.
#' @export
as_cell_program_data <- function(pseudobulk, sample_metadata = NULL,
                                 features = c("all", "variable")) {
  if (is.character(features) && length(features) > 1L) features <- features[1L]

  if (!is.list(pseudobulk) || is.null(names(pseudobulk)) ||
      any(names(pseudobulk) == "") || anyDuplicated(names(pseudobulk))) {
    .stopf("`pseudobulk` must be a named list with unique, non-empty names.")
  }
  for (ct in names(pseudobulk)) {
    m <- pseudobulk[[ct]]
    if (!is.matrix(m) || !is.numeric(m)) {
      .stopf("Cell type '%s': entry must be a numeric matrix.", ct)
    }
    if (is.null(rownames(m)) || anyDuplicated(rownames(m))) {
      .stopf("Cell type '%s': matrix must have unique row names (observation IDs).", ct)
    }
    if (is.null(colnames(m)) || anyDuplicated(colnames(m))) {
      .stopf("Cell type '%s': matrix must have unique column names (gene IDs).", ct)
    }
  }

  if (!is.null(sample_metadata)) {
    if (!is.data.frame(sample_metadata) || is.null(sample_metadata$observation_id)) {
      .stopf("`sample_metadata` must be a data.frame with an `observation_id` column.")
    }
    if (anyDuplicated(sample_metadata$observation_id)) {
      .stopf("`sample_metadata$observation_id` contains duplicates.")
    }
    all_ids <- unique(unlist(lapply(pseudobulk, rownames)))
    missing_meta <- setdiff(all_ids, sample_metadata$observation_id)
    if (length(missing_meta) > 0L) {
      .warnf("%d observation IDs have no metadata row (e.g. '%s').",
             length(missing_meta), missing_meta[1L])
    }
  }

  if (is.list(features)) {
    if (!identical(sort(names(features)), sort(names(pseudobulk)))) {
      .stopf("Named `features` list must have the same names as `pseudobulk`.")
    }
    for (ct in names(features)) {
      bad <- setdiff(features[[ct]], colnames(pseudobulk[[ct]]))
      if (length(bad) > 0L) {
        .stopf("Cell type '%s': %d requested features absent from the matrix (e.g. '%s').",
               ct, length(bad), bad[1L])
      }
    }
  } else {
    features <- match.arg(features, c("all", "variable"))
  }

  structure(
    list(
      matrices = pseudobulk,
      sample_metadata = sample_metadata,
      features = features,
      cell_types = names(pseudobulk)
    ),
    class = "cell_program_data"
  )
}

#' Validate a cell_program_data object
#'
#' Hard errors for structural problems; warnings for statistical concerns
#' (small N, near-zero-variance genes, sparse coverage).
#'
#' @param x A `cell_program_data` object.
#' @param min_obs Warn when a cell type has fewer than this many observations.
#' @return `x`, invisibly.
#' @export
validate_cell_program_data <- function(x, min_obs = 20L) {
  if (!inherits(x, "cell_program_data")) {
    .stopf("`x` must be a cell_program_data object.")
  }
  for (ct in x$cell_types) {
    m <- x$matrices[[ct]]
    if (nrow(m) < min_obs) {
      .warnf("Cell type '%s': only %d observations (< %d); factor estimates may be unstable.",
             ct, nrow(m), min_obs)
    }
    gv <- apply(m, 2L, stats::var, na.rm = TRUE)
    if (any(gv < .Machine$double.eps, na.rm = TRUE)) {
      .warnf("Cell type '%s': %d near-zero-variance genes.", ct,
             sum(gv < .Machine$double.eps, na.rm = TRUE))
    }
    frac_na <- mean(is.na(m))
    if (frac_na > 0.2) {
      .warnf("Cell type '%s': %.0f%% missing entries.", ct, 100 * frac_na)
    }
  }
  invisible(x)
}

#' @export
print.cell_program_data <- function(x, ...) {
  cat("cell_program_data\n")
  cat("  cell types:", length(x$cell_types), "\n")
  for (ct in x$cell_types) {
    m <- x$matrices[[ct]]
    cat(sprintf("    %-12s %4d obs x %5d genes (%.0f%% missing)\n",
                ct, nrow(m), ncol(m), 100 * mean(is.na(m))))
  }
  invisible(x)
}

#' @export
summary.cell_program_data <- function(object, ...) {
  data.frame(
    cell_type = object$cell_types,
    n_obs = vapply(object$matrices, nrow, integer(1L)),
    n_genes = vapply(object$matrices, ncol, integer(1L)),
    frac_missing = vapply(object$matrices, function(m) mean(is.na(m)), numeric(1L)),
    row.names = NULL
  )
}
