## cand_subspace.R -- subspace_sum joint-subspace caller (reviewed prototype on winner-consistent prep).
## A coordinated program is a donor direction held simultaneously by several cell types' score subspaces.
## Operator: A = Dh (sum_c P_c) Dh on the union of donors, P_c = projector of cell type c's (rank-first,
## site-residualized) score subspace with zero rows for donors c lacks, Dh = 1/sqrt(#cell types observing
## each donor). Top-J generalized eigen-directions (M v = lambda D v) are tested against within-strata row
## permutations of each cell type's basis (row permutation preserves U'U = I, so each cell type keeps its
## subspace but cross-cell-type donor alignment is broken). Per-direction p = (1 + #{lam_perm >= lam_obs}) /
## (n_perm + 1) with cummax step-down. Membership per direction by contrib[c] = ||U_c' v||^2 / ||v||^2 over
## c's observed donors, gated on its permutation distribution:
##   gate = "ct95": per-cell-type 95% quantile (reviewed prototype; uncorrected across cell types)
##   gate = "maxT": 95% quantile of the per-permutation max over cell types (FWER-ish within direction)
## Member program = score column with max |cor| to v. Pinned per review: J=15, n_perm=199, alpha=0.05,
## min_shared=20. Prep identical to n_agg_noanc (rank-first + ambient drop 0.3) so the arm isolates the
## caller statistic. Needs candidates.R (.basis, .mcca_setup, .mcca_op), cand_edge.R (.prep_rank_first),
## ambient.R (ambient_scores), methods.R (.perm_rows, .empty_clusters).

.subspace_spec <- function(Ub, st, J) {
  A <- .mcca_op(Ub, st)
  ev <- eigen(A, symmetric = TRUE)
  JJ <- min(J, length(ev$values))
  V <- ev$vectors[, seq_len(JJ), drop = FALSE] * st$Dh        # generalized eigenvectors in donor space
  contrib <- vapply(names(Ub), function(c) {
    o <- st$obs[[c]]; v <- V[o, , drop = FALSE]
    colSums(crossprod(Ub[[c]], v)^2) / pmax(colSums(v^2), 1e-12)
  }, numeric(JJ))
  colnames(contrib) <- names(Ub)
  list(values = ev$values[seq_len(JJ)], V = V, contrib = contrib)   # contrib: JJ x n_ct
}

cand_subspace <- function(fit, strata, seed = 1L, J = 15L, n_perm = 199L, alpha = 0.05,
                          min_shared = 20L, amb_thr = 0.3, gate = c("maxT", "ct95"), resid_site = TRUE) {
  set.seed(seed)
  gate <- match.arg(gate)
  Z <- .prep_rank_first(fit, strata, resid_site)
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {   # same ambient drop as n_agg_noanc
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  if (length(Z) < 2L) return(.empty_clusters())
  Ub <- lapply(Z, .basis)
  st <- .mcca_setup(Ub)
  ob <- .subspace_spec(Ub, st, J); J <- length(ob$values)
  nl <- matrix(NA_real_, n_perm, J); nc <- array(NA_real_, c(n_perm, J, length(Ub)))
  for (b in seq_len(n_perm)) {
    Up <- lapply(names(Ub), function(c) .perm_rows(Ub[[c]], strata[rownames(Ub[[c]])])); names(Up) <- names(Ub)
    s <- .subspace_spec(Up, st, J); nl[b, ] <- s$values; nc[b, , ] <- s$contrib
  }
  p <- cummax(vapply(seq_len(J), function(j) (1 + sum(nl[, j] >= ob$values[j])) / (n_perm + 1), 1))
  nobs <- vapply(Z, nrow, 1L)
  out <- list(); id <- 0L
  for (j in seq_len(J)) {
    if (p[j] >= alpha) break
    thr <- if (gate == "ct95") apply(nc[, j, , drop = FALSE], 2L, stats::quantile, 0.95, na.rm = TRUE)
           else stats::quantile(apply(nc[, j, , drop = FALSE], 1L, max), 0.95, na.rm = TRUE)
    mem <- which(ob$contrib[j, ] > thr & nobs[colnames(ob$contrib)] >= min_shared)
    if (length(mem) < 2L) next
    v <- stats::setNames(ob$V[, j], st$ids)
    ctm <- colnames(ob$contrib)[mem]
    prog <- vapply(ctm, function(ct) {
      S <- Z[[ct]]
      r <- suppressWarnings(abs(stats::cor(S, v[rownames(S)]))); r[!is.finite(r)] <- 0
      colnames(S)[which.max(r)]
    }, "")
    id <- id + 1L
    out[[id]] <- data.frame(cluster = id, cell_type = ctm, program_id = unname(prog),
                            stat = ob$values[j], fwer_p = p[j], stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else .empty_clusters()
}
