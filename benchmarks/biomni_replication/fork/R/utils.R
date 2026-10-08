## utils.R — internal helpers

`%||%` <- function(a, b) if (is.null(a)) b else a

.rms <- function(x) sqrt(mean(x^2))

.program_ids <- function(ct, k) sprintf("%s__F%03d", ct, seq_len(k))

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
