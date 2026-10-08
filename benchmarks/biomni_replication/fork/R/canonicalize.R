## canonicalize.R — resolve scale and sign indeterminacy
##
## For each factor k: scale so that sample scores have unit RMS, absorbing the
## scale into gene loadings; flip sign so the largest-|loading| gene is
## positive (deterministic tie-break by gene order). Reconstruction Z W' is
## exactly preserved.

#' Canonicalize factor scale and sign
#'
#' @param fit A `cell_program_fit` object.
#' @param tol Reconstruction-invariance tolerance (relative, Frobenius).
#' @return The fit with canonicalized scores/loadings (in place).
#' @export
canonicalize_programs <- function(fit, tol = 1e-6) {
  if (!inherits(fit, "cell_program_fit")) .stopf("`fit` must be a cell_program_fit.")
  if (isTRUE(fit$canonicalization$applied)) {
    .warnf("Programs already canonicalized; returning unchanged.")
    return(fit)
  }
  scales <- list()
  for (ct in names(fit$fits)) {
    Z <- fit$scores[[ct]]; W <- fit$loadings[[ct]]
    K <- ncol(Z)
    if (is.null(K) || K == 0L) { scales[[ct]] <- numeric(0); next }
    Z0 <- Z; W0 <- W
    s_k <- rep(1, K); sign_k <- rep(1L, K)
    for (k in seq_len(K)) {
      s <- .rms(Z[, k])
      if (s < .Machine$double.eps) next
      Z[, k] <- Z[, k] / s
      W[, k] <- W[, k] * s
      s_k[k] <- s
      g <- which.max(abs(W[, k]))
      if (length(g) > 0L && W[g, k] < 0) {
        Z[, k] <- -Z[, k]; W[, k] <- -W[, k]
        sign_k[k] <- -1L
      }
    }
    ## reconstruction invariance check
    R0 <- Z0 %*% t(W0)
    num <- max(abs(Z %*% t(W) - R0))
    den <- max(1e-12, sqrt(sum(R0^2)))
    if (num / den > tol) {
      .stopf("Canonicalization broke reconstruction in cell type '%s' (rel err %.2e).",
             ct, num / den)
    }
    fit$scores[[ct]] <- Z
    fit$loadings[[ct]] <- W
    scales[[ct]] <- data.frame(program_id = colnames(Z), scale = s_k,
                               sign_flip = sign_k, stringsAsFactors = FALSE)
  }
  fit$canonicalization <- list(applied = TRUE, scales = scales)
  fit
}
