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
#' @export
program_metadata <- function(fit) .program_table(fit)

#' Program score matrix across all cell types
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param format `"wide"` (observation x program, NA where a cell type was not
#'   observed; never imputed), `"long"`, or `"list"` (per-cell-type matrices).
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
