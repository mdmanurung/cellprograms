# MR-1: project fitted programs onto individual cells.
#
# Single-cell query data is put into the fit's feature space by replaying the
# exact preprocessing recorded at fit time (per-gene OLS residualization
# coefficients, then the centering/scale values actually used), and the
# canonical program activities are recovered by applying the same RMS/sign
# multipliers canonicalize_programs() applied to the donor scores. The
# underlying inference is flashier.utils::project_factor_activity() on the
# stored native fit; every check here fails closed.

#' Project fitted programs onto individual cells
#'
#' Estimates per-cell program activities for a `cell_program_fit`: each cell's
#' expression vector is adjusted with the fit's recorded preprocessing
#' (covariate residualization coefficients, then the stored centering and
#' scale values), and activities are inferred against the frozen gene
#' loadings with `flashier.utils::project_factor_activity()`. Raw activities
#' are canonicalized with the same RMS-scale and sign multipliers as the
#' donor scores, so projected activities are on the canonical program scale
#' and comparable across cell types and cohorts.
#'
#' The fit must be canonicalized (`canonicalize_programs()`) and must carry a
#' flashier.utils manifest (the default backend) with declared measurement
#' units (`feature_scale`); legacy-backend and fixed-covariate fits are
#' refused because their adjustment cannot be replayed. When the fit
#' residualized covariates for a cell type, `donors` and `sample_metadata`
#' are required so the adjustment can be rebuilt for the query cells; the
#' rebuilt design must reproduce the fit's design columns exactly (fail
#' closed on cohort changes that add or drop covariate levels).
#'
#' @param fit A canonicalized `cell_program_fit` (from
#'   `fit_celltype_programs()` + `canonicalize_programs()`).
#' @param expression Dense numeric matrix with cells in rows and genes in
#'   columns (row and column names required). Sparse matrices are not
#'   supported in v1. Genes absent from a cell type's fit panel are handled
#'   by `flashier.utils` intersection matching against `min_coverage`.
#' @param cell_types Named character vector mapping every row of
#'   `expression` to its cell type; names must be exactly the row names of
#'   `expression`. Cell types not in the fit (and fit cell types with no
#'   cells here) yield `NA` activities.
#' @param feature_scale Declared measurement units of `expression` (e.g.
#'   "logcounts"). Defaults to the fit's recorded `feature_scale`, which must
#'   not be "unknown"; a query whose declared units differ from the fit's is
#'   refused.
#' @param donors Optional named character vector mapping every cell to its
#'   donor ID (`sample_metadata$observation_id` values). Required whenever
#'   any projected cell type was fit with residualized covariates.
#' @param sample_metadata Optional data.frame with an `observation_id` column
#'   plus the covariate columns used at fit time. Required with `donors`
#'   when replaying residualization.
#' @param min_coverage Minimum fraction of a cell type's fit genes that must
#'   be present in `expression` (passed to the projection oracle).
#' @param maxiter Maximum projection iterations per cell type.
#' @param query_var_type Residual variance model for the query data: 0 =
#'   constant (default for single-cell profiles), 1 = per-cell, 2 = per-gene.
#' @param dense_budget_bytes Memory budget for dense intermediates (passed to
#'   the projection oracle).
#' @param keep_native When TRUE, include the native flash fit object for each
#'   projected cell type in `per_type[[ct]]$native`.
#'
#' @return A list with components:
#'   `activity` (cells x programs matrix; canonicalized activities, `NA` for
#'   cells whose type contributed no programs), `per_type` (per cell type:
#'   `activity`, `reconstruction` metrics, and optionally `native`),
#'   `coverage` (per cell type: programs, cells, shared/reference genes,
#'   coverage), `analysis_metadata`, and `provenance`.
#'
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable",
#'                              n_variable_genes = 150, feature_scale = "simulated")
#' fit <- canonicalize_programs(fit)
#' ct1 <- fit$cell_types[1]
#' # project the training pseudobulk back through the fit (donors as "cells")
#' cells <- setNames(rep(ct1, nrow(x$pseudobulk[[ct1]])), rownames(x$pseudobulk[[ct1]]))
#' proj <- project_programs_to_cells(fit, x$pseudobulk[[ct1]], cell_types = cells)
#' cor(proj$activity[, 1], fit$scores[[ct1]][, 1])   # high by construction
#' @export
project_programs_to_cells <- function(fit, expression, cell_types,
                                      feature_scale = NULL, donors = NULL,
                                      sample_metadata = NULL,
                                      min_coverage = 0.9, maxiter = 200L,
                                      query_var_type = 0L,
                                      dense_budget_bytes = 256 * 1024^2,
                                      keep_native = FALSE) {
  stopifnot(inherits(fit, "cell_program_fit"))
  if (is.null(fit$scores) || is.null(fit$canonicalization)) {
    cli::cli_abort("{.arg fit} must be canonicalized first: run canonicalize_programs(fit).")
  }
  if (!is.matrix(expression) || !is.numeric(expression)) {
    cli::cli_abort("{.arg expression} must be a dense numeric matrix with cells in rows and genes in columns; sparse inputs are not supported in v1.")
  }
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    cli::cli_abort("{.arg expression} must carry cell row names and gene column names.")
  }
  if (!is.character(cell_types) || is.null(names(cell_types)) ||
      !setequal(names(cell_types), rownames(expression))) {
    cli::cli_abort("{.arg cell_types} must be a named character vector mapping every row of {.arg expression} to its cell type.")
  }
  if (!is.null(donors) &&
      (!is.character(donors) || is.null(names(donors)) ||
       !all(names(cell_types) %in% names(donors)))) {
    cli::cli_abort("{.arg donors} must be a named character vector mapping every cell to its donor ID.")
  }
  if (!is.null(sample_metadata) &&
      !"observation_id" %in% names(sample_metadata)) {
    cli::cli_abort("{.arg sample_metadata} must contain an {.val observation_id} column.")
  }
  fs <- if (!is.null(feature_scale)) feature_scale else fit$preprocessing$feature_scale
  if (is.null(fs) || !is.character(fs) || length(fs) != 1L ||
      identical(fs, "unknown")) {
    cli::cli_abort(c(
      "The fit's measurement units are undeclared ({.val unknown}), so query compatibility cannot be checked.",
      i = "Refit with {.code feature_scale = \"logcounts\"} (or pass {.arg feature_scale} here)."
    ))
  }
  if (!requireNamespace("flashier.utils", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg flashier.utils} is required for projection.")
  }

  cells <- names(cell_types)
  unmodeled <- setdiff(unique(cell_types), fit$cell_types)
  if (length(unmodeled)) {
    cli::cli_warn("Cell type(s) {.val {unmodeled}} are not in the fit; their program activities are NA.")
  }
  prog_ids <- unlist(lapply(fit$cell_types, function(ct) {
    K <- ncol(fit$scores[[ct]])
    if (K > 0L) paste0(ct, "_", seq_len(K)) else character(0)
  }), use.names = FALSE)
  activity <- matrix(NA_real_, nrow = length(cells), ncol = length(prog_ids),
                     dimnames = list(cells, prog_ids))

  per_type <- list()
  per_type_meta <- list()
  cov_rows <- list()
  for (ct in fit$cell_types) {
    K <- ncol(fit$scores[[ct]])
    cells_ct <- cells[cell_types == ct]
    row <- data.frame(cell_type = ct, programs = K, cells = length(cells_ct),
                      reference_genes = length(fit$fits[[ct]]$genes),
                      shared_genes = NA_integer_, coverage = NA_real_)
    cov_rows[[length(cov_rows) + 1L]] <- row
    if (K == 0L) next          # skipped/degenerate cell type: no programs
    if (length(cells_ct) == 0L) next
    f <- fit$fits[[ct]]$flash
    if (!is.null(f$fixed_L) || is.null(attr(f, "flashbridge_metadata"))) {
      cli::cli_abort("Cell type {.val {ct}}: the stored fit has no flashier.utils manifest (legacy or fixed-covariate backend) and its adjustment cannot be replayed; refit with the default backend to project.")
    }
    replay <- fit$fits[[ct]]$preproc_replay
    if (is.null(replay)) {
      cli::cli_abort("Cell type {.val {ct}}: the fit predates preprocessing replay records; refit to project.")
    }
    genes_ct <- fit$fits[[ct]]$genes
    shared <- intersect(genes_ct, colnames(expression))
    if (!length(shared)) {
      cli::cli_abort("Cell type {.val {ct}}: none of the fit's {length(genes_ct)} genes are present in {.arg expression}.")
    }
    # Replay the fit's adjustment on the query, in the fit's order:
    # per-gene OLS residualization, then centering, then scaling.
    if (!is.null(replay$residualize_coefficients)) {
      if (is.null(donors) || is.null(sample_metadata)) {
        cli::cli_abort("Cell type {.val {ct}} was fit with residualized covariates: {.arg donors} and {.arg sample_metadata} are required to replay the adjustment on query cells.")
      }
      X <- .replay_design(donors[cells_ct], sample_metadata,
                          .spec_covariates(fit$preprocessing$covariates),
                          replay$design_columns, ct)
      beta <- replay$residualize_coefficients
      # R2/F2 (covariate audit): aliased columns of a rank-deficient design
      # carry NA coefficients (lm.fit pivot); the in-sample fitted values used
      # the pivoted solution, i.e. NA <=> 0.
      beta[is.na(beta)] <- 0
      Q <- as.matrix(expression[cells_ct, shared, drop = FALSE]) -
        X %*% beta[replay$design_columns, shared, drop = FALSE]
    } else {
      Q <- as.matrix(expression[cells_ct, shared, drop = FALSE])
    }
    idx <- match(shared, genes_ct)
    Q <- sweep(Q, 2, replay$center[idx], "-")
    Q <- sweep(Q, 2, replay$scale[idx], "/")

    pr <- flashier.utils::project_factor_activity(
      fit = f, query = Q, sample_side = "rows", feature_scale = fs,
      preprocessing = list(scope = "cellprograms_projection_replay"),
      factors = NULL, ebnm_fn = NULL, maxiter = maxiter,
      feature_match = "intersection", min_coverage = min_coverage,
      var_type = query_var_type, dense_budget_bytes = dense_budget_bytes)

    A <- sweep(sweep(as.matrix(pr$representation$activity), 2,
                     fit$canonicalization[[ct]]$scale_multiplier, "/"),
               2, fit$canonicalization[[ct]]$sign_multiplier, "*")
    colnames(A) <- paste0(ct, "_", seq_len(K))
    activity[cells_ct, colnames(A)] <- A
    row$shared_genes <- length(shared)
    row$coverage <- pr$analysis_metadata$coverage
    cov_rows[[length(cov_rows)]] <- row
    per_type[[ct]] <- list(
      activity = A,
      reconstruction = pr$reconstruction,
      native = if (keep_native) pr$fit else NULL
    )
    per_type_meta[[ct]] <- pr$analysis_metadata
  }

  self_version <- tryCatch(
    as.character(utils::packageVersion("cellprograms")), error = function(e) NULL)
  list(
    activity = activity,
    per_type = per_type,
    coverage = do.call(rbind, cov_rows),
    analysis_metadata = list(
      feature_scale = fs,
      query_var_type = query_var_type,
      min_coverage = min_coverage,
      maxiter = maxiter,
      unmodeled_cell_types = unmodeled,
      per_type = per_type_meta
    ),
    provenance = list(
      timestamp = Sys.time(),
      cellprograms_version = self_version,
      flashier_utils_version = tryCatch(
        as.character(utils::packageVersion("flashier.utils")), error = function(e) NULL),
      feature_scale = fs
    )
  )
}

# Rebuild the covariate design for the query donors of one cell type and
# check it against the design recorded at fit time. Fail closed: any column
# difference (a covariate level lost or gained among query donors) aborts,
# because silently projecting through a different adjustment would corrupt
# the activities.
.replay_design <- function(donors, sample_metadata, covariates,
                           design_columns, ct) {
  absent <- setdiff(unique(donors), sample_metadata$observation_id)
  if (length(absent)) {
    cli::cli_abort("{length(absent)} query donor(s) missing from {.arg sample_metadata}, e.g. {.val {absent[1]}}.")
  }
  f <- if (inherits(covariates, "formula")) covariates
       else stats::reformulate(covariates)
  vars <- all.vars(f)
  m <- sample_metadata[match(donors, sample_metadata$observation_id),
                       vars, drop = FALSE]
  if (anyNA(m)) {
    cli::cli_abort("Missing covariate values for query donors of cell type {.val {ct}}; impute or drop them first.")
  }
  X <- stats::model.matrix(f, m)
  if (!identical(colnames(X), design_columns)) {
    cli::cli_abort(c(
      "Replayed covariate design for cell type {.val {ct}} does not match the design recorded at fit time.",
      i = "A covariate level present at fit time is absent among the query donors (or a new level appears): the adjustment cannot be replayed. Align the query cohort or refit."
    ))
  }
  X
}
