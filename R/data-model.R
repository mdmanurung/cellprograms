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
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' x
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
    # CAND-03: reject ALL non-finite entries (NA/NaN/Inf), not just NA — Inf
    # previously passed this check and crashed validate/fit with an opaque
    # "missing value where TRUE/FALSE needed" error.
    if (any(!is.finite(m))) {
      n_bad <- sum(!is.finite(m))
      cli::cli_abort(
        "Element {.val {ct}} contains {n_bad} non-finite value(s) (NA/NaN/Inf); \\
         pseudobulk matrices must be complete and finite."
      )
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

  if (!is.null(feature_metadata)) {
    if (!is.data.frame(feature_metadata)) {
      cli::cli_abort("{.arg feature_metadata} must be a data.frame keyed by gene ID.")
    }
    genes_all <- unique(unlist(lapply(pseudobulk, colnames), use.names = FALSE))
    keyed <- if (!is.null(feature_metadata$gene_id)) {
      feature_metadata$gene_id
    } else if (!is.null(feature_metadata$gene)) {
      feature_metadata$gene
    } else {
      rownames(feature_metadata)
    }
    if (length(intersect(keyed, genes_all)) == 0L) {
      cli::cli_warn("{.arg feature_metadata} shares no gene IDs with the pseudobulk matrices; check its key (rownames or a {.val gene}/{.val gene_id} column).")
    }
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
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' validate_cell_program_data(x)
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
    # n = 1: var() is NA for every gene; na.rm keeps the small-N warning above
    # as the signal instead of crashing `if (n_zero > 0)` on an NA.
    n_zero <- sum(gene_var < .Machine$double.eps, na.rm = TRUE)
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
#' @param x A `cell_program_data` object.
#' @param ... Ignored.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' print(x)
#' @keywords internal
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
#' @param object A `cell_program_data` object.
#' @param ... Ignored.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' summary(x)
#' @keywords internal
#' @export
summary.cell_program_data <- function(object, ...) {
  print(object)
  invisible(NULL)
}
