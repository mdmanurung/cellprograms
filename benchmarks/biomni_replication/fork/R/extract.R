## extract.R — program scores, loadings, summaries, top genes

#' Program score matrix across all cell types
#'
#' @param fit A `cell_program_fit` object.
#' @param format `"wide"` (observation x program matrix aligned by
#'   observation_id, NA where a cell type was not observed), `"long"`
#'   (data.frame), or `"list"` (per-cell-type matrices).
#' @return Scores in the requested format. Missing values are always explicit
#'   NA — never silently imputed.
#' @export
program_scores <- function(fit, format = c("wide", "long", "list")) {
  format <- match.arg(format)
  if (format == "list") return(fit$scores)

  ids <- unique(unlist(fit$input$observation_ids))
  progs <- fit$programs$program_id
  Z <- matrix(NA_real_, nrow = length(ids), ncol = length(progs),
              dimnames = list(ids, progs))
  for (ct in names(fit$scores)) {
    S <- fit$scores[[ct]]
    if (ncol(S) == 0L) next
    Z[rownames(S), colnames(S)] <- S
  }
  if (format == "wide") return(Z)

  long <- do.call(rbind, lapply(names(fit$scores), function(ct) {
    S <- fit$scores[[ct]]
    if (ncol(S) == 0L) return(NULL)
    data.frame(
      observation_id = rep(rownames(S), times = ncol(S)),
      cell_type = ct,
      program_id = rep(colnames(S), each = nrow(S)),
      score = as.vector(S),
      stringsAsFactors = FALSE
    )
  }))
  long
}

#' Program loading matrices
#'
#' @param fit A `cell_program_fit` object.
#' @param format `"list"` (per-cell-type gene x program matrices) or `"long"`.
#' @export
program_loadings <- function(fit, format = c("list", "long")) {
  format <- match.arg(format)
  if (format == "list") return(fit$loadings)
  do.call(rbind, lapply(names(fit$loadings), function(ct) {
    W <- fit$loadings[[ct]]
    if (ncol(W) == 0L) return(NULL)
    data.frame(
      gene = rep(rownames(W), times = ncol(W)),
      cell_type = ct,
      program_id = rep(colnames(W), each = nrow(W)),
      loading = as.vector(W),
      stringsAsFactors = FALSE
    )
  }))
}

#' Program metadata table
#' @param fit A `cell_program_fit` object.
#' @export
program_metadata <- function(fit) fit$programs

#' Summary of programs per cell type
#' @param fit A `cell_program_fit` object.
#' @export
program_summary <- function(fit) {
  per_ct <- lapply(names(fit$fits), function(ct) {
    f <- fit$fits[[ct]]
    Y <- fit$input$matrices[[ct]]
    K <- f$n_factors
    R2 <- NA_real_
    if (K > 0L) {
      Yhat <- fit$scores[[ct]] %*% t(fit$loadings[[ct]])
      resid <- Y - Yhat
      R2 <- 1 - sum(resid^2, na.rm = TRUE) /
        sum((Y - rowMeans(Y, na.rm = TRUE))^2, na.rm = TRUE)
    }
    data.frame(
      cell_type = ct,
      n_obs = nrow(Y),
      n_genes = ncol(Y),
      n_programs = K,
      R2 = R2,
      elbo = f$elbo %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, per_ct)
}

#' Top genes of a program
#'
#' @param fit A `cell_program_fit` object.
#' @param program Program ID, e.g. `"Mono__F003"`.
#' @param n Number of genes.
#' @param by Rank by `"abs_loading"` or `"lfsr"` (local false sign rate).
#' @export
top_genes <- function(fit, program, n = 20L, by = c("abs_loading", "lfsr")) {
  by <- match.arg(by)
  ct <- fit$programs$cell_type[match(program, fit$programs$program_id)]
  if (is.na(ct)) .stopf("Unknown program '%s'.", program)
  W <- fit$loadings[[ct]][, program]
  out <- data.frame(
    gene = names(W),
    loading = as.numeric(W),
    stringsAsFactors = FALSE
  )
  lf <- fit$posteriors$lfsr[[ct]]
  if (!is.null(lf) && program %in% colnames(lf)) {
    out$lfsr <- as.numeric(lf[, program])
  }
  ord <- if (by == "abs_loading") order(-abs(out$loading)) else order(out$lfsr, -abs(out$loading))
  head(out[ord, ], n)
}
