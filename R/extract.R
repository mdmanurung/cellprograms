# cellprograms: extraction helpers (scores, loadings, summaries, top genes).
# All take a canonicalized `cell_program_fit`; program ids are `<cell_type>_<k>`.

.program_table <- function(fit) {
  do.call(rbind, lapply(names(fit$loadings), function(ct) {
    p <- colnames(fit$loadings[[ct]])
    if (!length(p)) return(NULL)
    data.frame(program_id = p, cell_type = ct, stringsAsFactors = FALSE)
  }))
}

#' Program metadata table (one row per program)
#' @param fit A canonicalized `cell_program_fit`.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' program_metadata(fit)
#' @export
program_metadata <- function(fit) .program_table(fit)

#' Program score matrix across all cell types
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param format `"wide"` (observation x program, NA where a cell type was not
#'   observed; never imputed), `"long"`, or `"list"` (per-cell-type matrices).
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' Z <- program_scores(fit)   # observations x programs, NA where unobserved
#' dim(Z)
#' head(program_scores(fit, "long"))
#' @export
program_scores <- function(fit, format = c("wide", "long", "list")) {
  format <- match.arg(format)
  if (format == "list") return(fit$scores)
  if (format == "long") {
    return(do.call(rbind, lapply(names(fit$scores), function(ct) {
      S <- fit$scores[[ct]]
      if (!ncol(S)) return(NULL)
      data.frame(observation_id = rep(rownames(S), ncol(S)), cell_type = ct,
                 program_id = rep(colnames(S), each = nrow(S)),
                 score = as.vector(S), stringsAsFactors = FALSE)
    })))
  }
  progs <- .program_table(fit)$program_id
  Z <- matrix(NA_real_, length(fit$observation_ids), length(progs),
              dimnames = list(fit$observation_ids, progs))
  for (S in fit$scores) if (ncol(S)) Z[rownames(S), colnames(S)] <- S
  Z
}

#' Program loading matrices
#' @param fit A canonicalized `cell_program_fit`.
#' @param format `"list"` (gene x program per cell type) or `"long"`.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' W <- program_loadings(fit)
#' sapply(W, dim)
#' @export
program_loadings <- function(fit, format = c("list", "long")) {
  format <- match.arg(format)
  if (format == "list") return(fit$loadings)
  do.call(rbind, lapply(names(fit$loadings), function(ct) {
    W <- fit$loadings[[ct]]
    if (!ncol(W)) return(NULL)
    data.frame(gene = rep(rownames(W), ncol(W)), cell_type = ct,
               program_id = rep(colnames(W), each = nrow(W)),
               loading = as.vector(W), stringsAsFactors = FALSE)
  }))
}

#' Summary of programs per cell type (n, K, R2 of the reconstruction, ELBO)
#' @param fit A canonicalized `cell_program_fit`.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' program_summary(fit)
#' @export
program_summary <- function(fit) {
  do.call(rbind, lapply(names(fit$fits), function(ct) {
    Y <- fit$fits[[ct]]$Y_used
    K <- ncol(fit$scores[[ct]])
    R2 <- NA_real_
    if (K > 0L) {
      res <- Y - fit$scores[[ct]] %*% t(fit$loadings[[ct]])
      R2 <- 1 - sum(res^2) / sum(sweep(Y, 2, colMeans(Y))^2)
    }
    data.frame(cell_type = ct, n_obs = nrow(Y), n_genes = ncol(Y),
               n_programs = K, R2 = R2,
               elbo = fit$fits[[ct]]$flash$elbo %||% NA_real_,
               stringsAsFactors = FALSE)
  }))
}

#' Top genes of a program
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param program Program ID, e.g. `"Mono_3"`.
#' @param n Number of genes.
#' @param by Rank by `"abs_loading"` or `"lfsr"` (flashier local false sign
#'   rate; needs the fit's gene-side lfsr).
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' top_genes(fit, "B_1", n = 5)
#' @export
top_genes <- function(fit, program, n = 20L, by = c("abs_loading", "lfsr")) {
  by <- match.arg(by)
  tab <- .program_table(fit)
  ct <- tab$cell_type[match(program, tab$program_id)]
  if (is.na(ct)) .stopf("Unknown program '%s'.", program)
  W <- fit$loadings[[ct]][, program]
  out <- data.frame(gene = names(W), loading = as.numeric(W), stringsAsFactors = FALSE)
  lf <- fit$fits[[ct]]$flash$F_lfsr
  if (!is.null(lf)) out$lfsr <- as.numeric(lf[, match(program, colnames(fit$loadings[[ct]]))])
  if (by == "lfsr" && is.null(lf)) .stopf("No lfsr stored for cell type '%s'.", ct)
  ord <- if (by == "abs_loading") order(-abs(out$loading)) else order(out$lfsr, -abs(out$loading))
  utils::head(out[ord, ], n)
}

#' Refit program loadings on a common (union) gene panel
#'
#' When each cell type selected its own features (e.g. `features = "variable"`),
#' the per-cell-type loading matrices live on different gene sets and cannot be
#' compared across cell types. This helper takes the union of the fitted gene
#' panels and, for each cell type, regresses the union-panel expression on that
#' cell type's program scores by least squares. The fits and scores are not
#' modified: this is route B of ROADMAP item C5 (the `genes_after` arm of the
#' model-selection benchmark showed that fitting on the union panel from the
#' start, route A, does not improve recovery). Genes a cell type did not select
#' get a least-squares loading on its programs instead of being missing; for
#' the genes a cell type did fit, the refit loadings are close to (but not
#' identical to) the flashier posterior means, which are shrunk.
#'
#' The fit's preprocessing is replayed on the panel before the regression:
#' covariate residualization (residualize mode) or removal of the fitted
#' covariate part (fixed mode), then the recorded centering/scaling.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param genes Optional character vector restricting the panel; default is
#'   the union of the per-cell-type fitted genes. Genes absent from a cell
#'   type's pseudobulk get `NA` loadings there, as do genes that are degenerate
#'   (zero-variance after centering) in that cell type when scaling.
#' @return A list with `loadings`, a named list of gene x program matrices on
#'   the common panel (one per cell type, 0 columns for skipped cell types),
#'   and `genes`, the panel used.
#'
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3,
#'                                                    features = "variable",
#'                                                    n_variable_genes = 150))
#' up <- union_panel_loadings(fit)
#' sapply(up$loadings, nrow)   # same gene rows for every cell type
#' up$genes[1:5]
#' @export
union_panel_loadings <- function(fit, genes = NULL) {
  if (!inherits(fit, "cell_program_fit")) {
    .stopf("`fit` must be a cell_program_fit (run fit_celltype_programs() first).")
  }
  panel <- if (is.null(genes)) {
    unique(unlist(lapply(fit$fits, function(f) f$genes), use.names = FALSE))
  } else genes
  if (!is.character(panel) || !length(panel)) {
    .stopf("`genes` must be NULL or a non-empty character vector.")
  }
  pp <- fit$preprocessing
  loadings <- stats::setNames(lapply(names(fit$fits), function(ct) {
    S <- fit$scores[[ct]]
    K <- if (is.null(S)) 0L else ncol(S)
    G <- matrix(NA_real_, length(panel), K,
                dimnames = list(panel, if (is.null(S)) NULL else colnames(S)))
    if (!K) return(G)
    avail <- intersect(panel, colnames(fit$pseudobulk[[ct]]))
    Y <- as.matrix(fit$pseudobulk[[ct]][rownames(S), avail, drop = FALSE])
    # Replay the fit's covariate handling: residualize (residualize mode) or
    # subtract the OLS covariate part under the fit's own fixed design (fixed
    # mode), so the regression on the scores sees the same signal the fit did.
    if (!is.null(pp$covariates) && pp$covariate_mode == "residualize") {
      X <- .cov_design(rownames(S), fit$sample_metadata, pp$covariates, ct)
      if (!is.null(X)) Y <- stats::lm.fit(X, Y)$residuals
    } else {
      Xc <- fit$fits[[ct]]$flash$fixed_L
      if (!is.null(Xc)) Y <- Y - Xc %*% solve(crossprod(Xc), crossprod(Xc, Y))
    }
    if (isTRUE(pp$centered)) Y <- scale(Y, center = TRUE, scale = isTRUE(pp$scaled))
    else if (isTRUE(pp$scaled)) Y <- scale(Y, center = FALSE, scale = TRUE)
    G[avail, ] <- t(stats::lm.fit(S, Y)$coefficients)   # genes x programs
    G[is.nan(G)] <- NA_real_   # degenerate panel genes under scaling
    G
  }), names(fit$fits))
  list(loadings = loadings, genes = panel)
}
