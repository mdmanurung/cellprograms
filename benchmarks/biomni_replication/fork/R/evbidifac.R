## evbidifac.R — Level 3 candidate B: EV-BIDIFAC-style linked decomposition
##
## Independent implementation of the empirical variational Bayes linked
## matrix decomposition published as:
##   Lock EF. "Empirical Bayes Linked Matrix Decomposition."
##   Machine Learning 113:7451-7477 (2024). arXiv:2408.00237.
##
## The algorithm (manuscript Algorithm 1, Theorems 6-7) is implemented here
## from the published equations. The authors' reference scripts
## (github.com/lockEF/bidifac, no explicit license) are NOT shipped with this
## package; they are used only as a numerical test oracle during development.
##
## Model: X_ij = sum_k U_i^(k) V_j^(k)' + E_ij over modules k defined by
## row-set/column-set indicators; module structure is estimated by EVB
## singular-value shrinkage with per-block noise variances.

## ---------------------------------------------------------------------------
## Single-matrix EVB pieces (manuscript Theorems 6-7)
## ---------------------------------------------------------------------------

## Zero-crossing kappa of the EVB threshold equation, by dense grid search.
.evb_kappa <- function(M, N, lb = 0.5, ub = 4, n_grid = 5000L) {
  alpha <- M / N
  k <- seq(lb, ub, length.out = n_grid)
  f <- log(sqrt(alpha) * k + 1) / (sqrt(alpha) * k) +
       log(k / sqrt(alpha) + 1) / (k / sqrt(alpha)) - 1
  k[which.min(f^2)]
}

## Free-energy objective psi for noise-variance estimation (Theorem 7).
.evb_psi <- function(x, alpha, kappa) {
  xbar <- 1 + alpha + sqrt(alpha) * (kappa + 1 / kappa)
  res <- x - log(x)
  ind <- which(x >= xbar)
  if (length(ind) > 0L) {
    kf <- (1 / (2 * sqrt(alpha))) *
      ((x[ind] - (1 + alpha)) + sqrt(pmax((x[ind] - (1 + alpha))^2 - 4 * alpha, 0)))
    res[ind] <- x[ind] - log(x[ind]) + log(sqrt(alpha) * kf + 1) +
      alpha * log(kf / sqrt(alpha) + 1) - sqrt(alpha) * kf
  }
  res
}

## Estimate noise variance sigma^2 of a single matrix (Theorem 7).
.evb_sigma2 <- function(X, kappa = NULL, lb = 0.01, n_grid = 500L) {
  M <- max(dim(X)); N <- min(dim(X))
  alpha <- N / M
  if (is.null(kappa)) kappa <- .evb_kappa(M, N)
  d <- svd(X)$d
  d <- d[d > 1e-3]
  if (length(d) == 0L) return(mean(X^2))
  ub <- mean(X^2)
  if (!is.finite(ub) || ub <= 0) return(1)
  sig_vec <- seq(lb * ub, ub, length.out = n_grid)
  om <- vapply(sig_vec, function(s2) {
    sum(.evb_psi(d^2 / (M * s2), alpha, kappa))
  }, numeric(1L))
  sig_vec[which.min(om)]
}

## EVB singular-value shrinkage of a single matrix (Theorem 6).
## sig2 is the noise variance; returns the denoised low-rank estimate.
.evb_shrink <- function(X, sig2 = 1, kappa = NULL) {
  M <- nrow(X); N <- ncol(X)
  if (is.null(kappa)) kappa <- .evb_kappa(M, N)
  thr <- sqrt(sig2) * sqrt(M + N + sqrt(M * N) * (kappa + 1 / kappa))
  s <- svd(X)
  d <- s$d
  keep <- d > thr
  d_new <- numeric(length(d))
  if (any(keep)) {
    dk <- d[keep]
    d_new[keep] <- (dk / 2) * (1 - (M + N) * sig2 / dk^2 +
      sqrt(pmax((1 - (M + N) * sig2 / dk^2)^2 - 4 * M * N * sig2^2 / dk^4, 0)))
  }
  list(S = s$u %*% diag(d_new, nrow = length(d_new)) %*% t(s$v),
       rank = sum(d_new > 1e-10),
       d = d_new, u = s$u, v = s$v)
}

## Nuclear-norm soft threshold (warm start), threshold sqrt(M) + sqrt(N).
.nn_shrink <- function(X) {
  M <- nrow(X); N <- ncol(X)
  s <- svd(X)
  thr <- sqrt(M) + sqrt(N)
  d_new <- pmax(s$d - thr, 0)
  s$u %*% diag(d_new, nrow = length(d_new)) %*% t(s$v)
}

## ---------------------------------------------------------------------------
## Full decomposition over modules
## ---------------------------------------------------------------------------

#' EV-BIDIFAC-style linked decomposition
#'
#' @param X Matrix (rows x columns), possibly with NA.
#' @param col_sets List of column-index vectors (e.g., one per cell type).
#' @param row_sets List of row-index vectors (default: single set, all rows).
#' @param modules List of modules, each `list(rows = <row-set ids>,
#'   cols = <col-set ids>)`. Default: all non-empty row-set x col-set
#'   combinations (the full BIDIFAC+ module enumeration).
#' @param max_iter Maximum cyclic iterations (EVB phase).
#' @param warm_iter Maximum warm-start iterations (nuclear-norm phase). The
#'   warm phase must run to convergence — stopping early leaves structure in
#'   small modules and biases the EVB phase away from the global module.
#' @param conv_tol Convergence tolerance (relative change in fitted values).
#' @param verbose Print iteration progress.
#' @return A list with per-module reconstructed structures (original scale),
#'   per-block noise variances, and module metadata.
#' @export
evbidifac_decompose <- function(X, col_sets, row_sets = NULL,
                                modules = NULL, max_iter = 200L,
                                warm_iter = 500L, conv_tol = 1e-6,
                                verbose = FALSE) {
  if (is.null(row_sets)) row_sets <- list(seq_len(nrow(X)))
  I <- length(row_sets); J <- length(col_sets)

  if (is.null(modules)) {
    ## enumerate all non-empty row-set x col-set combinations
    modules <- list()
    for (rsub in .powerset(seq_len(I))) {
      for (csub in .powerset(seq_len(J))) {
        modules[[length(modules) + 1L]] <- list(rows = rsub, cols = csub)
      }
    }
  }
  K <- length(modules)

  ## module submatrix indices
  mod_rows <- lapply(modules, function(m) unlist(row_sets[m$rows]))
  mod_cols <- lapply(modules, function(m) unlist(col_sets[m$cols]))

  ## observed mask; impute NA with 0 initially
  isna <- is.na(X)
  Ximp <- X; Ximp[isna] <- 0

  ## per-block noise variance (on observed entries)
  sigma2 <- matrix(1, I, J)
  for (i in seq_len(I)) for (j in seq_len(J)) {
    blk <- X[row_sets[[i]], col_sets[[j]], drop = FALSE]
    if (any(!is.na(blk))) {
      sigma2[i, j] <- if (anyNA(blk)) mean(blk^2, na.rm = TRUE) else .evb_sigma2(blk)
    }
  }

  ## scale blocks to unit noise variance
  Xs <- Ximp
  for (i in seq_len(I)) for (j in seq_len(J)) {
    Xs[row_sets[[i]], col_sets[[j]]] <- Xs[row_sets[[i]], col_sets[[j]]] / sqrt(sigma2[i, j])
  }

  ## module structures on the scaled scale + running total fit
  S <- lapply(seq_len(K), function(k) matrix(0, length(mod_rows[[k]]), length(mod_cols[[k]])))
  Fitted <- matrix(0, nrow(X), ncol(X))
  Fitted[isna] <- 0

  cycle_modules <- function(shrink_fn, n_iter) {
    for (iter in seq_len(n_iter)) {
      delta <- 0
      for (k in seq_len(K)) {
        rr <- mod_rows[[k]]; cc <- mod_cols[[k]]
        R_sub <- Xs[rr, cc, drop = FALSE] - Fitted[rr, cc, drop = FALSE] + S[[k]]
        new <- shrink_fn(R_sub)
        delta <- delta + sum((new - S[[k]])^2)
        Fitted[rr, cc] <<- Fitted[rr, cc] + (new - S[[k]])
        S[[k]] <<- new
      }
      ## EM imputation of missing entries
      if (any(isna)) Xs[isna] <<- Fitted[isna]
      if (verbose) message(sprintf("  iter %d: delta = %.4g", iter, delta))
      if (delta < conv_tol * max(1, sum(Fitted^2))) break
    }
    iter
  }

  ## warm start (nuclear norm) then EVB phase
  n_warm <- cycle_modules(.nn_shrink, warm_iter)
  n_evb <- cycle_modules(function(R) .evb_shrink(R, sig2 = 1)$S, max_iter)

  ## NOTE: no post-hoc sigma2 update. The Theorem-7 per-block estimates are
  ## the paper's estimator and are accurate; a residual-based update must be
  ## computed on a consistent scale (Fitted is on the scaled scale), and the
  ## reference implementation (glob.bidi.o) does not update sigma2 either.

  ## unscale modules to the original scale and summarize
  mod_info <- data.frame(
    module_id = sprintf("M%02d", seq_len(K)),
    row_sets = vapply(modules, function(m) paste(m$rows, collapse = "+"), character(1L)),
    col_sets = vapply(modules, function(m) paste(m$cols, collapse = "+"), character(1L)),
    n_col_sets = vapply(modules, function(m) length(m$cols), integer(1L)),
    class = vapply(modules, function(m) {
      if (length(m$cols) == J) "global"
      else if (length(m$cols) == 1L) paste0("private:", m$cols)
      else paste0("partial:", paste(m$cols, collapse = "+"))
    }, character(1L)),
    rank = integer(K),
    frob = numeric(K),
    stringsAsFactors = FALSE
  )

  S_out <- vector("list", K)
  for (k in seq_len(K)) {
    ## unscale: multiply each block of module k by sigma_ij
    S_un <- S[[k]]
    rlab <- rep(modules[[k]]$rows, times = lengths(row_sets[modules[[k]]$rows]))
    clab <- rep(modules[[k]]$cols, times = lengths(col_sets[modules[[k]]$cols]))
    for (i in unique(rlab)) for (j in unique(clab)) {
      S_un[rlab == i, clab == j] <- S_un[rlab == i, clab == j] * sqrt(sigma2[i, j])
    }
    S_out[[k]] <- list(rows = mod_rows[[k]], cols = mod_cols[[k]], S = S_un)
    mod_info$rank[k] <- sum(svd(S_un)$d > 1e-8)
    mod_info$frob[k] <- sqrt(sum(S_un^2))
  }

  list(
    modules = mod_info,
    S = S_out,
    sigma2 = sigma2,
    n_iter = c(warm = n_warm, evb = n_evb),
    converged = n_evb < max_iter
  )
}

.powerset <- function(x) {
  n <- length(x)
  out <- list()
  for (m in seq_len(n)) {
    out <- c(out, utils::combn(x, m, simplify = FALSE))
  }
  out
}

## ---------------------------------------------------------------------------
## Package-facing wrapper: EV-BIDIFAC as an integration method
## ---------------------------------------------------------------------------

#' EV-BIDIFAC integration of cell-type programs (Level 3 candidate B)
#'
#' @param fit A `cell_program_fit`.
#' @param level `"scores"` (decompose Z_all) or `"expression"` (stacked
#'   cell-type expression matrices).
#' @param max_modules Optional cap on the module enumeration (default: all
#'   non-empty cell-type subsets).
#' @param seed Unused (algorithm is deterministic); accepted for interface
#'   consistency.
#' @param ... Passed to `evbidifac_decompose`.
#' @return An integration-style list with method `"ev_bidifac"`.
#' @export
fit_joint_evbidifac <- function(fit, level = c("scores", "expression"),
                                max_modules = NULL, seed = NULL, ...) {
  level <- match.arg(level)
  st <- .build_stacked(fit, level)
  X <- st$X; block_id <- st$block_id
  cts <- unique(block_id)
  col_sets <- lapply(cts, function(ct) which(block_id == ct))
  names(col_sets) <- cts

  modules <- NULL
  if (!is.null(max_modules)) {
    ps <- .powerset(seq_along(cts))
    modules <- lapply(utils::head(ps, max_modules), function(csub) list(rows = 1L, cols = csub))
  }

  dec <- evbidifac_decompose(X, col_sets = col_sets, modules = modules, ...)

  ## translate set-index class labels to cell-type names
  set_names <- names(col_sets)
  if (!is.null(set_names)) {
    remap <- function(x) {
      if (x == "global") return("global")
      parts <- strsplit(x, ":", fixed = TRUE)[[1L]]
      idx <- as.integer(strsplit(parts[2L], "+", fixed = TRUE)[[1L]])
      paste0(parts[1L], ":", paste(set_names[idx], collapse = "+"))
    }
    dec$modules$class <- vapply(dec$modules$class, remap, character(1L))
  }

  ## drop empty modules; build states + weights from module SVDs
  keep <- which(dec$modules$frob > 1e-8)
  states <- list(); weights <- list(); state_meta <- list()
  for (k in keep) {
    Sk <- dec$S[[k]]$S
    rownames(Sk) <- rownames(X)[dec$S[[k]]$rows]
    colnames(Sk) <- colnames(X)[dec$S[[k]]$cols]
    sv <- svd(Sk)
    r <- sum(sv$d > 1e-8)
    if (r == 0L) next
    for (rr in seq_len(r)) {
      sid <- sprintf("%s_r%d", dec$modules$module_id[k], rr)
      states[[sid]] <- sv$u[, rr] * sv$d[rr]
      weights[[sid]] <- sv$v[, rr]
      state_meta[[sid]] <- list(
        module = dec$modules$module_id[k],
        class = dec$modules$class[k],
        col_sets = dec$modules$col_sets[k],
        variance = sv$d[rr]^2
      )
    }
  }

  ## build states/weights matrices embedded in the full row/column spaces
  ## (module submatrices have different widths; zero-pad outside the module)
  S_mat <- matrix(NA_real_, nrow(X), length(states),
                  dimnames = list(rownames(X), names(states)))
  W_mat <- matrix(0, ncol(X), length(states),
                  dimnames = list(colnames(X), names(states)))
  for (k in keep) {
    Sk <- dec$S[[k]]$S
    sv <- svd(Sk); r <- sum(sv$d > 1e-8)
    if (r == 0L) next
    rids <- rownames(X)[dec$S[[k]]$rows]
    cids <- colnames(X)[dec$S[[k]]$cols]
    for (rr in seq_len(r)) {
      sid <- sprintf("%s_r%d", dec$modules$module_id[k], rr)
      S_mat[rids, sid] <- sv$u[, rr] * sv$d[rr]
      W_mat[cids, sid] <- sv$v[, rr]
    }
  }

  out <- list(
    method = "ev_bidifac",
    level = level,
    decomposition = dec,
    states = S_mat,
    state_meta = state_meta,
    block_id = block_id,
    block_scales = st$block_scales
  )
  if (level == "scores") {
    out$program_weights <- W_mat
  } else {
    gl <- list()
    for (ct in cts) {
      W <- W_mat[block_id == ct, , drop = FALSE]
      rownames(W) <- sub(paste0("^", ct, "__"), "", rownames(W))
      gl[[ct]] <- W
    }
    out$gene_loadings <- gl
  }
  out$summary <- data.frame(
    state = names(states),
    module = vapply(state_meta, `[[`, character(1L), "module"),
    class = vapply(state_meta, `[[`, character(1L), "class"),
    variance = vapply(state_meta, `[[`, numeric(1L), "variance"),
    stringsAsFactors = FALSE
  )
  out
}
