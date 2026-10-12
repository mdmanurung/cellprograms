## candidates.R -- candidate CMP callers. Every candidate: function(fit, strata, seed, ...) -> clusters df
## (cluster, cell_type, program_id, stat, fwer_p). Needs methods.R (.nodes, .perm_rows, .empty_clusters).

## ---- shared helpers ---------------------------------------------------------------------
## centred (optionally site-residualized) node scores per cell type, dropping empty cell types
.prep_scores <- function(fit, strata, resid_site = TRUE, min_donors = 20L) {
  Z <- fit$scores[vapply(fit$scores, function(S) !is.null(S) && ncol(S) > 0 && nrow(S) >= min_donors, TRUE)]
  lapply(Z, function(S) {
    f <- factor(strata[rownames(S)])
    X <- if (resid_site && nlevels(f) > 1) stats::model.matrix(~f) else matrix(1, nrow(S), 1)
    R <- S - X %*% qr.solve(X, S); dimnames(R) <- dimnames(S); R
  })
}

## orthonormal basis of the column space of a centred score matrix
.basis <- function(S, tol = 1e-8) {
  s <- svd(S); k <- sum(s$d > tol * s$d[1]); U <- s$u[, seq_len(k), drop = FALSE]; rownames(U) <- rownames(S); U
}

## MCCA-style operator: A = Dh (sum_c P_c) Dh on the union of donors, P_c = projector of cell type c
## (zero rows for donors c lacks); Dd = # cell types observing each donor.
.mcca_setup <- function(Ub) {
  ids <- sort(unique(unlist(lapply(Ub, rownames))))
  obs <- lapply(Ub, function(U) match(rownames(U), ids))
  Dd <- tabulate(unlist(obs), length(ids)); Dh <- 1 / sqrt(pmax(Dd, 1))
  list(ids = ids, obs = obs, Dh = Dh, n = length(ids))
}
.mcca_op <- function(Ub, st, drop = NULL) {
  A <- matrix(0, st$n, st$n)
  for (c in setdiff(names(Ub), drop)) { o <- st$obs[[c]]; A[o, o] <- A[o, o] + tcrossprod(Ub[[c]]) }
  st$Dh * A * rep(st$Dh, each = st$n)
}
.top_eig <- function(A, k = 1L) {
  if (k == 1L && nrow(A) > 60) {
    e <- tryCatch(RSpectra::eigs_sym(A, 1L, which = "LA"), error = function(e) NULL)
    if (!is.null(e)) return(list(values = e$values, vectors = e$vectors))
  }
  e <- eigen(A, symmetric = TRUE); list(values = e$values[seq_len(k)], vectors = e$vectors[, seq_len(k), drop = FALSE])
}

## ---- C1: leave-one-cell-type-out MCCA with node-removal step-down ------------------------
## global test: lambda_1 of the MCCA operator vs within-site row permutations of every cell type
## membership: node k of cell type c is tested against the consensus built WITHOUT c (v_{-c}), which is
## independent of c under H0(c) -> Fisher-z test, Holm across all (c, k). Members' nodes are removed and
## the procedure repeats (max J components).
cand_loo_mcca <- function(fit, strata, seed = 1L, n_perm = 199L, alpha = 0.05, alpha_mem = 0.05, J = 6L,
                          resid_site = TRUE, min_cts = 2L, mem_rule = "holm", rg = FALSE) {
  set.seed(seed)
  Z <- .prep_scores(fit, strata, resid_site)
  if (rg) Z <- lapply(Z, .rank_gauss)
  out <- list(); id <- 0L
  for (j in seq_len(J)) {
    Z <- Z[vapply(Z, ncol, 1L) > 0]
    if (length(Z) < min_cts) break
    Ub <- lapply(Z, .basis); st <- .mcca_setup(Ub)
    lam <- .top_eig(.mcca_op(Ub, st))$values[1]
    nl <- vapply(seq_len(n_perm), function(b) {
      Up <- lapply(names(Ub), function(c) .perm_rows(Ub[[c]], strata[rownames(Ub[[c]])])); names(Up) <- names(Ub)
      .top_eig(.mcca_op(Up, st))$values[1]
    }, 1)
    p <- (1 + sum(nl >= lam)) / (n_perm + 1)
    if (p >= alpha) break
    ## leave-one-out consensus and node tests
    tests <- do.call(rbind, lapply(names(Z), function(c) {
      v <- .top_eig(.mcca_op(Ub, st, drop = c))$vectors[, 1] * st$Dh
      names(v) <- st$ids
      S <- Z[[c]]; vv <- v[rownames(S)]
      r <- suppressWarnings(as.numeric(stats::cor(S, vv))); r[!is.finite(r)] <- 0
      n <- nrow(S)
      zst <- atanh(pmin(abs(r), 0.999999)) * sqrt(n - 3)
      data.frame(cell_type = c, program_id = colnames(S), r = r, z = zst, p = 2 * stats::pnorm(-zst), stringsAsFactors = FALSE)
    }))
    tests$padj <- stats::p.adjust(tests$p, if (mem_rule == "holm") "holm" else "BH")
    mem <- tests[tests$padj < alpha_mem, , drop = FALSE]
    if (length(unique(mem$cell_type)) < min_cts) break
    id <- id + 1L
    out[[id]] <- data.frame(cluster = id, cell_type = mem$cell_type, program_id = mem$program_id, stat = lam, fwer_p = p,
                            stringsAsFactors = FALSE)
    for (k in seq_len(nrow(mem))) Z[[mem$cell_type[k]]] <- Z[[mem$cell_type[k]]][, colnames(Z[[mem$cell_type[k]]]) != mem$program_id[k], drop = FALSE]
  }
  if (length(out)) do.call(rbind, out) else .empty_clusters()
}
