## integrate-pmd.R — optional stage-2 integration: multi-view penalized
## matrix decomposition (PMD) of per-cell-type program scores.
##
## DIALOGUE-style coupling (Jerby-Arnon & Regev, Nat Biotechnol 2022): each
## cell type's stage-1 EBMF program-score matrix is one "view"; find sparse
## unit-norm program-weight vectors a_1..a_C maximizing the sum of pairwise
## correlations of view scores across samples (Witten, Tibshirani & Hastie
## 2009 multi-view PMD / sparse CCA, SUMCOR criterion). L1 penalties are
## tuned by permutation; successive multicellular programs (MCPs) are
## extracted under deflation, each with an empirical p-value from an
## independent permutation batch.
##
## Missingness: sample x cell-type coverage is never imputed. Views keep
## their own observed samples; the objective sums PAIRWISE-COMPLETE
## correlations (each view pair evaluated on samples observed in both), and
## the alternating updates use available-case targets (each sample
## contributes only the views in which it is observed). Views are
## standardized on their own observed samples.
##
## The solver implements the Witten-Tibshirani alternating soft-threshold
## updates directly (no PMA dependency), allowing single-program views
## (K_c = 1) and deterministic SVD initialization.
##
## Design contract (per package principles):
##   - integration never replaces stage-1 programs;
##   - full traceability: MCP -> program weights -> gene loadings
##     (trace_state()).

## ---------------------------------------------------------------------------
## Solver
## ---------------------------------------------------------------------------

## Map shared relative sparsity rho in (0,1] to per-view L1 bounds.
## ||a_c||_1 ranges [1, sqrt(K_c)] for unit-norm a_c; smaller = sparser.
.pmd_sumabs <- function(X, rho) {
  Ks <- vapply(X, ncol, integer(1L))
  1 + rho * (sqrt(Ks) - 1)
}

## Soft-threshold operator.
.pmd_soft <- function(x, lam) sign(x) * pmax(abs(x) - lam, 0)

## Binary search for the smallest lambda achieving ||a||_1 <= bound
## (a = soft(g, lambda) / ||soft(g, lambda)||_2).
.pmd_lambda <- function(g, bound) {
  n2 <- sqrt(sum(g^2))
  if (n2 < .Machine$double.eps) return(0)
  if (sum(abs(g / n2)) <= bound) return(0)  # unconstrained already sparse enough
  lo <- 0; hi <- max(abs(g))
  for (i in seq_len(60L)) {
    mid <- (lo + hi) / 2
    a <- .pmd_soft(g, mid)
    n2 <- sqrt(sum(a^2))
    if (n2 < .Machine$double.eps) { lo <- mid; next }
    if (sum(abs(a / n2)) > bound) lo <- mid else hi <- mid
  }
  hi
}

## One multi-view PMD component at fixed per-view L1 bounds, with
## pairwise-complete (available-case) handling of missing sample x view
## coverage. X: named list of (n_c x K_c) matrices with sample-ID rownames.
## Returns per-view weights, per-view scores (on each view's own samples),
## the pairwise-complete SUMCOR objective (NA if degenerate), the matrix of
## pairwise sample counts, and the pairwise correlation matrix.
.pmd_solve <- function(X, sumabs, n_iter = 100L, tol = 1e-7, min_pair = 15L) {
  C <- length(X)
  ids <- lapply(X, rownames)
  Ks <- vapply(X, ncol, integer(1L))
  ## deterministic init: first right singular vector per view
  ws <- lapply(X, function(x) {
    v <- svd(x, nu = 0, nv = 1)$v[, 1L]
    v / sqrt(sum(v^2))
  })
  for (it in seq_len(n_iter)) {
    ws_old <- ws
    for (c in seq_len(C)) {
      ## available-case target: sum of other views' variates where observed
      t <- numeric(nrow(X[[c]]))
      for (d in setdiff(seq_len(C), c)) {
        ud <- as.numeric(X[[d]] %*% ws[[d]])
        m <- match(ids[[c]], ids[[d]])
        has <- !is.na(m)
        t[has] <- t[has] + ud[m[has]]
      }
      g <- as.numeric(crossprod(X[[c]], t))
      lam <- .pmd_lambda(g, sumabs[c])
      a <- .pmd_soft(g, lam)
      n2 <- sqrt(sum(a^2))
      ws[[c]] <- if (n2 > .Machine$double.eps) a / n2 else rep(0, Ks[c])
    }
    delta <- max(vapply(seq_len(C),
                        function(c) max(abs(ws[[c]] - ws_old[[c]])), numeric(1L)))
    if (delta < tol) break
  }
  U <- Map(function(x, w) as.numeric(x %*% w), X, ws)
  if (any(vapply(U, stats::sd, numeric(1L)) < 1e-10)) {
    return(list(ws = ws, U = U, objective = NA_real_))
  }
  ## pairwise-complete objective
  obj <- 0
  npair <- 0L
  pair_n <- matrix(NA_integer_, C, C, dimnames = list(names(X), names(X)))
  pair_r <- matrix(NA_real_, C, C, dimnames = list(names(X), names(X)))
  for (c in seq_len(C - 1L)) {
    for (d in (c + 1L):C) {
      m <- match(ids[[c]], ids[[d]])
      both <- !is.na(m)
      pair_n[c, d] <- pair_n[d, c] <- sum(both)
      if (sum(both) < min_pair) next
      r <- stats::cor(U[[c]][both], U[[d]][m[both]])
      pair_r[c, d] <- pair_r[d, c] <- r
      obj <- obj + r
      npair <- npair + 1L
    }
  }
  if (npair == 0L) return(list(ws = ws, U = U, objective = NA_real_))
  list(ws = ws, U = U, objective = obj, pair_n = pair_n, pair_r = pair_r)
}

## Permutation tuning of the sparsity grid (DIALOGUE-style): shuffle values
## within each view independently while keeping rownames fixed — because
## .pmd_solve aligns views by rowname, permuting rows together with their
## rownames would be a no-op (label<->value pairing survives). Shuffling
## values against fixed labels breaks cross-view coupling while preserving
## within-view structure and each view's sample set. Select rho by empirical
## p-value (ties broken by larger observed objective).
.pmd_permute_views <- function(X) {
  lapply(X, function(x) {
    x2 <- x[sample.int(nrow(x)), , drop = FALSE]
    rownames(x2) <- rownames(x)
    x2
  })
}
.pmd_tune <- function(X, rho_grid, n_perm, min_pair) {
  grid <- lapply(rho_grid, function(rho) .pmd_sumabs(X, rho))
  solve_obj <- function(Xx, sa) {
    o <- .pmd_solve(Xx, sa, min_pair = min_pair)$objective
    if (is.na(o)) -Inf else o
  }
  obs <- vapply(grid, function(sa) solve_obj(X, sa), numeric(1L))
  perm <- matrix(-Inf, n_perm, length(grid))
  for (b in seq_len(n_perm)) {
    Xp <- .pmd_permute_views(X)
    for (j in seq_along(grid)) {
      perm[b, j] <- solve_obj(Xp, grid[[j]])
    }
  }
  p <- (1 + colSums(sweep(perm, 2L, obs, `>=`))) / (1 + n_perm)
  best <- which(p == min(p))
  best <- best[which.max(obs[best])]
  list(sumabs = grid[[best]], rho = rho_grid[best],
       tuning_table = data.frame(rho = rho_grid, objective = obs,
                                 tuning_p = p, stringsAsFactors = FALSE))
}

## Independent permutation p-value at a fixed penalty (not used for tuning).
.pmd_perm_test <- function(X, sumabs, n_perm, min_pair) {
  obs <- .pmd_solve(X, sumabs, min_pair = min_pair)$objective
  if (is.na(obs)) return(list(objective = NA_real_, pvalue = NA_real_,
                              perm_null = numeric(0)))
  perm <- replicate(n_perm, {
    Xp <- .pmd_permute_views(X)
    .pmd_solve(Xp, sumabs, min_pair = min_pair)$objective
  })
  perm <- perm[!is.na(perm)]
  list(objective = obs,
       pvalue = (1 + sum(perm >= obs)) / (1 + length(perm)),
       perm_null = perm)
}

## Deflate each view along its canonical variate (PMD deflation), on the
## view's own samples.
.pmd_deflate <- function(X, ws) {
  Map(function(x, w) {
    u <- as.numeric(x %*% w)
    x - (tcrossprod(u) %*% x) / max(sum(u^2), .Machine$double.eps)
  }, X, ws)
}

## ---------------------------------------------------------------------------
## Main entry
## ---------------------------------------------------------------------------

#' Multi-view PMD integration of stage-1 programs
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param n_states Number of MCPs to extract (default min(5, min K_c)).
#' @param rho_grid Shared relative sparsity grid in (0, 1].
#' @param n_perm Permutations for tuning and for the per-MCP test (each).
#' @param min_pair Minimum pairwise-complete samples for a view pair to
#'   enter the objective.
#' @param standardize Standardize program scores per cell type (recommended:
#'   EBMF factor scores are on heterogeneous scales).
#' @param seed Random seed.
#' @param verbose Print progress.
#' @return Integration list (see integrate_programs).
.integrate_pmd <- function(fit, n_states = NULL,
                           rho_grid = c(0.15, 0.3, 0.5, 0.7, 1),
                           n_perm = 100L, min_pair = 15L,
                           standardize = TRUE, seed = NULL, verbose = 1L) {
  if (!is.null(seed)) set.seed(seed)

  scores <- program_scores(fit, "list")
  cts <- names(scores)[vapply(scores, function(S) ncol(S) > 0L, logical(1L))]
  if (length(cts) < 2L) .stopf("pmd integration requires >= 2 cell types with programs.")
  scores <- scores[cts]
  all_ids <- unique(unlist(lapply(scores, rownames)))

  ## drop zero-variance programs (on the view's own samples)
  for (ct in cts) {
    v <- apply(scores[[ct]], 2L, stats::sd)
    keep <- v >= .Machine$double.eps
    if (!all(keep)) {
      if (verbose > 0L) .warnf("Dropping %d zero-variance programs in '%s'.",
                               sum(!keep), ct)
      scores[[ct]] <- scores[[ct]][, keep, drop = FALSE]
    }
  }
  if (standardize) {
    scores <- lapply(scores, function(S) scale(S, center = TRUE, scale = TRUE))
  }

  ## pairwise coverage report
  C <- length(cts)
  pair_n0 <- matrix(NA_integer_, C, C, dimnames = list(cts, cts))
  for (c in seq_len(C - 1L)) {
    for (d in (c + 1L):C) {
      pair_n0[c, d] <- pair_n0[d, c] <-
        length(intersect(rownames(scores[[c]]), rownames(scores[[d]])))
    }
  }
  if (any(pair_n0 < min_pair, na.rm = TRUE)) {
    .warnf("Some view pairs have < %d common samples and are excluded from the objective.",
           min_pair)
  }

  K_max <- n_states %||% min(5L, min(vapply(scores, ncol, integer(1L))))
  K_max <- max(1L, K_max)

  Xd <- scores
  fnorm0 <- vapply(scores, function(x) sqrt(sum(x^2)), numeric(1L))
  mcps <- list()
  for (k in seq_len(K_max)) {
    ## drop deflation-exhausted views (e.g., single-program cell types whose
    ## only program was consumed by an earlier MCP)
    keep <- vapply(Xd, function(x) sqrt(sum(x^2)), numeric(1L)) >
      1e-6 * fnorm0[names(Xd)]
    if (sum(keep) < 2L) {
      if (verbose > 0L) {
        message(sprintf("MCP%02d: fewer than 2 non-exhausted views; stopping.", k))
      }
      break
    }
    if (!all(keep)) {
      if (verbose > 0L) {
        message(sprintf("MCP%02d: dropping exhausted view(s): %s",
                        k, paste(names(Xd)[!keep], collapse = ", ")))
      }
      Xd <- Xd[keep]
    }
    tuned <- .pmd_tune(Xd, rho_grid, n_perm, min_pair)
    sol <- .pmd_solve(Xd, tuned$sumabs, min_pair = min_pair)
    test <- .pmd_perm_test(Xd, tuned$sumabs, n_perm, min_pair)
    if (is.na(test$objective)) {
      if (verbose > 0L) {
        message(sprintf("MCP%02d: degenerate variate after deflation; stopping.", k))
      }
      break
    }

    ## per-view weights named by program ID
    ws <- Map(function(w, x) stats::setNames(w, colnames(x)), sol$ws, Xd)
    ## sign canonicalization: largest-|weight| program positive
    allw <- unlist(ws)
    flip <- if (allw[which.max(abs(allw))] < 0) -1 else 1
    ws <- lapply(ws, `*`, flip)

    ## per-view MCP scores for ALL samples, then NA-aware average
    active_cts <- names(ws)
    U_full <- Map(function(S, w) as.numeric(S %*% w), scores[active_cts], ws)
    U_mat <- matrix(NA_real_, length(all_ids), length(cts),
                    dimnames = list(all_ids, cts))
    for (ct in active_cts) {
      U_mat[rownames(scores[[ct]]), ct] <- U_full[[ct]]
    }
    state <- rowMeans(U_mat, na.rm = TRUE)

    ## participation diagnostics. SUMCOR forces unit-norm weights in every
    ## active view, so L1 weight share cannot localize an MCP; instead use
    ## each view's contribution to the objective (sum of pairwise cors with
    ## all other views).
    pair_r <- sol$pair_r
    diag(pair_r) <- 0
    contribution <- rowSums(pair_r, na.rm = TRUE)
    l1 <- vapply(ws, function(w) sum(abs(w)), numeric(1L))
    participation <- l1 / sum(l1)

    ## contribution-weighted state: downweights views that do not covary with
    ## the rest (private programs riding along under unit-norm weights)
    cw <- rep(0, length(cts)); names(cw) <- cts
    cw[active_cts] <- pmax(contribution[active_cts], 0)
    cw[is.na(cw)] <- 0
    state_weighted <- apply(U_mat, 1L, function(r) {
      ok <- !is.na(r)
      if (sum(cw[ok]) <= 0) return(mean(r[ok]))
      stats::weighted.mean(r[ok], cw[ok])
    })

    mcps[[k]] <- list(
      state_id = sprintf("MCP%02d", k),
      weights = ws,
      active_celltypes = active_cts,
      state = stats::setNames(state, all_ids),
      state_weighted = stats::setNames(state_weighted, all_ids),
      state_by_celltype = U_mat,
      pair_cors = pair_r,
      pair_n = sol$pair_n,
      view_contribution = contribution,
      objective = test$objective,
      pvalue = test$pvalue,
      rho = tuned$rho,
      tuning_table = tuned$tuning_table,
      participation = participation
    )
    if (verbose > 0L) {
      message(sprintf(paste0("MCP%02d: objective %.3f, perm p = %.3f, ",
                             "rho = %.2f, active CTs: %d"),
                      k, test$objective, test$pvalue, tuned$rho,
                      length(active_cts)))
    }
    Xd <- .pmd_deflate(Xd, sol$ws)
  }

  if (length(mcps) == 0L) {
    return(list(method = "pmd", states = NULL, program_weights = NULL,
                mcps = list(), cell_types = cts,
                summary = data.frame(method = "pmd"),
                note = "no non-degenerate MCPs found"))
  }

  ## package-contract outputs
  state_ids <- vapply(mcps, `[[`, character(1L), "state_id")
  S <- vapply(mcps, `[[`, numeric(length(all_ids)), "state")
  colnames(S) <- state_ids
  rownames(S) <- all_ids
  Sw <- vapply(mcps, `[[`, numeric(length(all_ids)), "state_weighted")
  colnames(Sw) <- state_ids
  rownames(Sw) <- all_ids
  prog_ids <- unlist(lapply(scores, colnames), use.names = FALSE)
  pw <- matrix(0, length(prog_ids), length(mcps),
               dimnames = list(prog_ids, state_ids))
  for (k in seq_along(mcps)) {
    for (ct in names(mcps[[k]]$weights)) {
      w <- mcps[[k]]$weights[[ct]]
      pw[names(w), k] <- w
    }
  }

  list(
    method = "pmd",
    states = S,
    states_weighted = Sw,
    states_by_celltype = stats::setNames(lapply(cts, function(ct) {
      m <- vapply(mcps, function(mcp) mcp$state_by_celltype[, ct],
                  numeric(length(all_ids)))
      rownames(m) <- all_ids
      colnames(m) <- state_ids
      m
    }), cts),
    program_weights = pw,
    mcps = mcps,
    cell_types = cts,
    pairwise_samples = pair_n0,
    standardize = standardize,
    rho_grid = rho_grid,
    n_perm = n_perm,
    min_pair = min_pair,
    summary = data.frame(
      method = "pmd",
      state = state_ids,
      objective = vapply(mcps, `[[`, numeric(1L), "objective"),
      perm_pvalue = vapply(mcps, `[[`, numeric(1L), "pvalue"),
      rho = vapply(mcps, `[[`, numeric(1L), "rho"),
      n_active_celltypes = vapply(mcps, function(m) length(m$active_celltypes),
                                  integer(1L)),
      active_celltypes = vapply(mcps, function(m) paste(m$active_celltypes, collapse = "+"),
                                character(1L)),
      stringsAsFactors = FALSE
    )
  )
}
