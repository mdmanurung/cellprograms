`%||%` <- function(a, b) if (is.null(a)) b else a

## utils.R — internal helpers

.stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

.warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

## Orthonormal basis for the column space of a matrix (rank-revealing via SVD)
.orthobasis <- function(W, tol = sqrt(.Machine$double.eps)) {
  if (ncol(W) == 0L) return(matrix(0, nrow(W), 0L))
  s <- svd(W)
  r <- sum(s$d > tol * max(s$d))
  if (r == 0L) return(matrix(0, nrow(W), 0L))
  s$u[, seq_len(r), drop = FALSE]
}

## Hungarian-style assignment via greedy matching on a similarity matrix.
## Greedy is adequate for the small factor counts involved (k <= max_factors).
.greedy_match <- function(sim) {
  sim <- abs(sim)
  # D-3b: NA/NaN in the similarity matrix made max(S) NA and crashed the loop
  # with "missing value where TRUE/FALSE needed". A non-finite similarity
  # carries no matching information: demote to -Inf (never selected) and say so.
  bad <- !is.finite(sim)
  if (any(bad)) {
    cli::cli_warn("{sum(bad)} non-finite similarity value(s) treated as unmatched.")
    sim[bad] <- -Inf
  }
  matches <- data.frame(a = integer(0), b = integer(0), sim = numeric(0))
  if (nrow(sim) == 0L || ncol(sim) == 0L) return(matches)
  ## track ORIGINAL row/col indices: submatrix coordinates shift after removal
  ri <- seq_len(nrow(sim)); ci <- seq_len(ncol(sim))
  S <- sim
  while (length(S) > 0 && max(S) > -Inf) {
    idx <- which(S == max(S), arr.ind = TRUE)[1L, ]
    matches <- rbind(matches, data.frame(a = ri[idx[1L]], b = ci[idx[2L]],
                                         sim = S[idx[1L], idx[2L]]))
    keep_r <- seq_along(ri) != idx[1L]; keep_c <- seq_along(ci) != idx[2L]
    S <- S[keep_r, keep_c, drop = FALSE]
    ri <- ri[keep_r]; ci <- ci[keep_c]
    if (length(S) == 0L) break
  }
  matches
}

# Fail-closed guard for fit consumers: fit$scores / fit$loadings per-type
# entries are populated by canonicalize_programs() only. On a raw fit these
# are empty lists, and indexing them crashes opaquely ("argument is of length
# zero" via ncol(NULL)) or returns silently-empty tables. The codebase already
# guards three entry points this way (coordination.R, cell-projection.R,
# plots.R); this helper makes the message uniform. The condition is the actual
# broken state (no populated scores/loadings) rather than the canonicalization
# marker alone, so minimal hand-constructed fits with valid matrices pass.
.require_canonicalized <- function(fit, fn) {
  has_scores <- !is.null(fit$scores) && length(fit$scores) > 0L
  has_loadings <- !is.null(fit$loadings) && length(fit$loadings) > 0L
  if (is.null(fit$canonicalization) && (!has_scores || !has_loadings)) {
    .stopf("%s requires a canonicalized fit; run canonicalize_programs() first.", fn)
  }
}
