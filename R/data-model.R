# cellprograms: data model
#
# The canonical input is a named list of observations x genes matrices
# (one per cell type) plus a sample metadata data.frame keyed by
# observation_id. Missing cell-type observations are allowed and are
# never silently imputed.

#' Create a cell program data object
#'
#' @param pseudobulk Named list of numeric matrices (observations x genes),
#'   one per cell type. Row names must be observation IDs.
#' @param sample_metadata A data.frame containing an `observation_id` column.
#'   Outcome variables may be present but are never used for factorization.
#' @param feature_metadata Optional data.frame keyed by gene ID.
#'
#' @return An object of class `cell_program_data`.
#'
#' @export
as_cell_program_data <- function(pseudobulk, sample_metadata = NULL,
                                 feature_metadata = NULL) {
  if (!is.list(pseudobulk) || length(pseudobulk) == 0L) {
    cli::cli_abort("{.arg pseudobulk} must be a non-empty named list of matrices.")
  }
  if (is.null(names(pseudobulk)) || any(!nzchar(names(pseudobulk)))) {
    cli::cli_abort("Every element of {.arg pseudobulk} must have a non-empty cell-type name.")
  }
  if (anyDuplicated(names(pseudobulk))) {
    cli::cli_abort("Duplicated cell-type names in {.arg pseudobulk}.")
  }

  for (ct in names(pseudobulk)) {
    m <- pseudobulk[[ct]]
    if (!is.matrix(m) || !is.numeric(m)) {
      cli::cli_abort("Element {.val {ct}} of {.arg pseudobulk} must be a numeric matrix.")
    }
    if (any(is.na(m))) {
      cli::cli_abort("Element {.val {ct}} contains NA values; pseudobulk matrices must be complete.")
    }
    if (is.null(rownames(m)) || any(!nzchar(rownames(m)))) {
      cli::cli_abort("Element {.val {ct}} must have non-empty row names (observation IDs).")
    }
    if (is.null(colnames(m)) || any(!nzchar(colnames(m)))) {
      cli::cli_abort("Element {.val {ct}} must have non-empty column names (gene IDs).")
    }
    if (anyDuplicated(rownames(m))) {
      cli::cli_abort("Duplicated observation IDs in cell type {.val {ct}}.")
    }
    if (anyDuplicated(colnames(m))) {
      cli::cli_abort("Duplicated gene IDs in cell type {.val {ct}}.")
    }
  }

  if (is.null(sample_metadata)) {
    sample_metadata <- data.frame(
      observation_id = unique(unlist(lapply(pseudobulk, rownames)))
    )
  }
  if (!is.data.frame(sample_metadata)) {
    cli::cli_abort("{.arg sample_metadata} must be a data.frame.")
  }
  if (!"observation_id" %in% names(sample_metadata)) {
    cli::cli_abort("{.arg sample_metadata} must contain an {.val observation_id} column.")
  }
  if (anyDuplicated(sample_metadata$observation_id)) {
    cli::cli_abort("Duplicated {.val observation_id} values in {.arg sample_metadata}.")
  }

  all_ids <- unique(unlist(lapply(pseudobulk, rownames)))
  unknown <- setdiff(all_ids, sample_metadata$observation_id)
  if (length(unknown) > 0L) {
    cli::cli_abort(
      "{length(unknown)} observation ID(s) present in pseudobulk matrices but
       missing from {.arg sample_metadata}: {.val {utils::head(unknown, 5)}}"
    )
  }

  # Gene consistency check: overlapping cell types should share gene IDs,
  # but different cell types may use different gene sets (per plan Section 8).
  # We record, not enforce, the intersection.
  gene_sets <- lapply(pseudobulk, colnames)
  shared_genes <- Reduce(intersect, gene_sets)

  obj <- list(
    pseudobulk = pseudobulk,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    cell_types = names(pseudobulk),
    observation_ids = all_ids,
    shared_genes = shared_genes,
    dims = lapply(pseudobulk, dim)
  )
  class(obj) <- "cell_program_data"
  obj
}

#' Validate a cell program data object
#'
#' Runs the full data-contract check and returns the object invisibly,
#' emitting warnings for soft-failure conditions (small N, near-zero
#' variance genes, extreme sparsity).
#'
#' @param x A `cell_program_data` object.
#'
#' @export
validate_cell_program_data <- function(x) {
  stopifnot(inherits(x, "cell_program_data"))
  for (ct in x$cell_types) {
    n <- x$dims[[ct]][1]
    if (n < 10) {
      cli::cli_warn(
        "Cell type {.val {ct}} has only {n} observations; EBMF is unreliable at very small N."
      )
    }
    m <- x$pseudobulk[[ct]]
    gene_var <- apply(m, 2, stats::var)
    n_zero <- sum(gene_var < .Machine$double.eps)
    if (n_zero > 0) {
      cli::cli_warn(
        "Cell type {.val {ct}}: {n_zero} gene(s) have near-zero variance."
      )
    }
    if (mean(m == 0) > 0.5) {
      cli::cli_warn(
        "Cell type {.val {ct}}: more than 50% of entries are exactly zero;
         this does not look like transformed pseudobulk expression."
      )
    }
  }
  invisible(x)
}

#' Print a cell program data object
#' @export
print.cell_program_data <- function(x, ...) {
  cli::cli_text("{.strong cell_program_data}")
  cli::cli_text("Cell types: {length(x$cell_types)} {.val {x$cell_types}}")
  for (ct in x$cell_types) {
    d <- x$dims[[ct]]
    cli::cli_text("  {ct}: {d[1]} observations x {d[2]} genes")
  }
  cli::cli_text("Shared genes across all cell types: {length(x$shared_genes)}")
  invisible(x)
}

#' Summarize a cell program data object
#' @export
summary.cell_program_data <- function(object, ...) {
  print(object)
  invisible(NULL)
}
