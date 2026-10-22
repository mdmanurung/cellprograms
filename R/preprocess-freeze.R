# MR-0: frozen preprocessing specifications.
#
# A frozen spec records the exact feature panel and adjustment decisions of
# an analysis once, so that refits (new cohorts, added donors, re-runs) use
# byte-identical preprocessing instead of silently re-deriving it. The spec
# is content-addressed (sha256 digest) and validated fail-closed by
# fit_celltype_programs(preprocessing_spec = ...).

#' Freeze a preprocessing specification
#'
#' Resolves the feature panel and adjustment decisions once and returns a
#' content-addressed specification (sha256 digest) that
#' `fit_celltype_programs(preprocessing_spec = ...)` replays exactly on any
#' dataset carrying the same genes. Intended workflow: freeze the spec on the
#' discovery cohort, store it alongside the fit (it is plain RDS-serializable
#' data), and pass it to every refit so the frozen panel can never silently
#' drift.
#'
#' @param x A `cell_program_data` object (the discovery cohort).
#' @param features Feature selection mode, passed to the same selection used
#'   by `fit_celltype_programs()`: "all", "variable", or a named list of gene
#'   ID vectors per cell type. The resolved per-cell-type gene vectors are
#'   stored in the spec.
#' @param n_variable_genes Number of most-variable genes when
#'   `features = "variable"`.
#' @param exclude_genes Optional character vector of gene IDs removed from
#'   every cell type's feature universe before feature selection.
#' @param exclude_artifacts Convenience preset: when TRUE, union
#'   `artifact_genes()` over the pseudobulk gene universe into
#'   `exclude_genes`. The spec stores the resolved union.
#' @param covariates Optional sample-level covariates: a character vector of
#'   `sample_metadata` column names or a one-sided formula (stored as its
#'   deparsed string so the spec serializes and digests stably).
#' @param covariate_mode How `covariates` are handled in refits:
#'   "residualize" (default) or "fixed".
#' @param center Center genes within each cell type (default TRUE).
#' @param scale Unit-variance scale genes (default FALSE; not recommended).
#' @param feature_scale Declared measurement units (default "unknown").
#'
#' @return A list of class `cell_program_preproc_spec` with components
#'   `spec_version`, `feature_sets` (named list of gene ID vectors per cell
#'   type), `exclude_genes` (resolved, sorted, unique), `exclude_artifacts`,
#'   `covariates` (character vector or deparsed formula string),
#'   `covariate_mode`, `center`, `scale`, `gene_universe` (sorted unique gene
#'   names of the discovery cohort), and `digest` (sha256 over the spec
#'   contents). The digest is what refits record in their preprocessing
#'   provenance.
#'
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 150)
#' fit <- fit_celltype_programs(x, max_factors = 3, preprocessing_spec = spec)
#' fit$preprocessing$frozen_spec == spec$digest   # TRUE
#' @export
freeze_preprocessing <- function(x, features = "all", n_variable_genes = 2000,
                                 exclude_genes = NULL, exclude_artifacts = FALSE,
                                 covariates = NULL,
                                 covariate_mode = c("residualize", "fixed"),
                                 center = TRUE, scale = FALSE,
                                 feature_scale = "unknown") {
  stopifnot(inherits(x, "cell_program_data"))
  covariate_mode <- match.arg(covariate_mode)
  stopifnot(is.logical(center), length(center) == 1L, !is.na(center),
            is.logical(scale), length(scale) == 1L, !is.na(scale),
            is.numeric(n_variable_genes), length(n_variable_genes) == 1L,
            n_variable_genes > 0)
  if (!is.null(exclude_genes) &&
      (!is.character(exclude_genes) || anyNA(exclude_genes))) {
    cli::cli_abort("{.arg exclude_genes} must be NULL or a character vector of gene IDs.")
  }
  universe <- unique(unlist(lapply(x$pseudobulk, colnames), use.names = FALSE))
  if (isTRUE(exclude_artifacts)) {
    exclude_genes <- union(exclude_genes, artifact_genes(universe))
  }
  if (!is.null(covariates)) {
    .validate_covariate_formula(covariates)
    vars <- if (inherits(covariates, "formula")) all.vars(covariates) else covariates
    absent <- setdiff(vars, names(x$sample_metadata))
    if (length(absent)) {
      cli::cli_abort("{.arg covariates} not in {.field sample_metadata}: {.val {absent}}.")
    }
  }
  feature_sets <- .select_features(x, features, n_variable_genes, exclude_genes)
  spec <- list(
    spec_version = "1.1.0",
    feature_sets = feature_sets,
    exclude_genes = if (is.null(exclude_genes)) NULL else sort(unique(exclude_genes)),
    exclude_artifacts = isTRUE(exclude_artifacts),
    covariates = if (inherits(covariates, "formula")) {
      paste(deparse(covariates), collapse = "")
    } else covariates,
    covariate_mode = covariate_mode,
    center = center,
    scale = scale,
    feature_scale = .validate_feature_scale(feature_scale),
    gene_universe = sort(universe)
  )
  class(spec) <- c("cell_program_preproc_spec", "list")
  spec$digest <- digest::digest(spec, algo = "sha256")
  spec
}

.validate_feature_scale <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    cli::cli_abort("{.arg feature_scale} must be a nonempty declared unit string.")
  }
  x
}

# Version 1.0 specs remain readable only when their original digest is intact.
# Their units are unknown; never infer units from numerical appearance.
.validate_preproc_spec <- function(spec) {
  if (!inherits(spec, "cell_program_preproc_spec")) {
    cli::cli_abort("{.arg preprocessing_spec} must come from freeze_preprocessing().")
  }
  if (!is.character(spec$spec_version) || length(spec$spec_version) != 1L ||
      is.na(spec$spec_version) || !spec$spec_version %in% c("1.0.0", "1.1.0") ||
      !is.character(spec$digest) || length(spec$digest) != 1L || is.na(spec$digest)) {
    cli::cli_abort("Unsupported preprocessing spec schema or missing digest; refreeze from the discovery data.")
  }
  original <- spec$digest
  spec$digest <- NULL
  if (!identical(original, digest::digest(spec, algo = "sha256"))) {
    cli::cli_abort("Frozen preprocessing spec digest mismatch: contents changed. Refreeze deliberately; do not reuse a stale hash.")
  }
  required <- c("feature_sets", "covariate_mode", "center", "scale", "gene_universe")
  if (!all(required %in% names(spec)) || !is.list(spec$feature_sets) ||
      is.null(names(spec$feature_sets)) || anyDuplicated(names(spec$feature_sets)) ||
      !length(spec$feature_sets) ||
      any(!vapply(spec$feature_sets, function(g) is.character(g) && length(g) > 0L &&
        !anyNA(g) && !anyDuplicated(g) && all(nzchar(g)), logical(1))) ||
      !spec$covariate_mode %in% c("residualize", "fixed") ||
      !is.logical(spec$center) || length(spec$center) != 1L || is.na(spec$center) ||
      !is.logical(spec$scale) || length(spec$scale) != 1L || is.na(spec$scale)) {
    cli::cli_abort("Malformed frozen preprocessing spec; refreeze from discovery data.")
  }
  if (spec$spec_version == "1.1.0") .validate_feature_scale(spec$feature_scale)
  invisible(TRUE)
}

# Rehydrate the covariates stored in a frozen spec: a deparsed one-sided
# formula string becomes a formula again; anything else passes through.
.spec_covariates <- function(cov) {
  if (is.character(cov) && length(cov) == 1L && grepl("^\\s*~", cov)) {
    return(stats::as.formula(cov))
  }
  cov
}
