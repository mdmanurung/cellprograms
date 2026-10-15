## cand_spec.R -- C2: node-level spectral (T-matrix) CMP caller with LOO membership + node-removal step-down.
## T_ij = atanh(r_ij) sqrt(n_ij - 3) for node pairs in different cell types (shared donors >= min_shared),
## 0 within cell type. Global stat = lambda_1(T) vs within-site row permutations of each cell type.
## Membership: for each cell type c, the leading eigenvector of T without c gives node weights u; the
## consensus v_{-c}(donor) = sum_i u_i z_i(donor) / sqrt(sum observed u_i^2); c's nodes tested vs v_{-c}
## (Fisher z, Holm/BH across all nodes). Member nodes are removed and the procedure repeats.
.std_cols <- function(S) { S <- scale(S); S[!is.finite(S)] <- 0; S }
.ct_pairs <- function(Z, min_shared = 20L) {
  cts <- names(Z); pairs <- list()
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    ids <- intersect(rownames(Z[[cts[i]]]), rownames(Z[[cts[j]]]))
    if (length(ids) >= min_shared) pairs[[length(pairs) + 1L]] <- list(a = cts[i], b = cts[j], ids = ids)
  }
  pairs
}
.node_T <- function(Z, pairs) {
  nodes <- unlist(lapply(names(Z), function(c) paste(c, colnames(Z[[c]]), sep = "|")))
  T <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (p in pairs) {
    r <- suppressWarnings(stats::cor(Z[[p$a]][p$ids, , drop = FALSE], Z[[p$b]][p$ids, , drop = FALSE])); r[!is.finite(r)] <- 0
    t <- atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(length(p$ids) - 3)
    ia <- paste(p$a, colnames(Z[[p$a]]), sep = "|"); ib <- paste(p$b, colnames(Z[[p$b]]), sep = "|")
    T[ia, ib] <- t; T[ib, ia] <- t(t)
  }
  T
}
## consensus donor score from all cell types except `c`, using node weights u (named "ct|prog")
.loo_consensus <- function(Z, u, c, d) {
  num <- numeric(length(d)); den <- numeric(length(d))
  for (o in setdiff(names(Z), c)) {
    w <- u[paste(o, colnames(Z[[o]]), sep = "|")]; w[is.na(w)] <- 0
    if (!any(w != 0)) next
    sc <- as.numeric(Z[[o]] %*% w); names(sc) <- rownames(Z[[o]])
    hit <- d %in% names(sc); num[hit] <- num[hit] + sc[d[hit]]; den[hit] <- den[hit] + sum(w^2)
  }
  ifelse(den > 0, num / sqrt(den), NA)
}
.node_tests <- function(S, vv) {
  ok <- is.finite(vv); n <- sum(ok)
  r <- if (n > 5) suppressWarnings(as.numeric(stats::cor(S[ok, , drop = FALSE], vv[ok]))) else rep(0, ncol(S))
  r[!is.finite(r)] <- 0
  zst <- atanh(pmin(abs(r), 0.999999)) * sqrt(max(n - 3, 1))
  data.frame(program_id = colnames(S), r = r, z = zst, p = 2 * stats::pnorm(-zst), n = n, stringsAsFactors = FALSE)
}

cand_node_spec <- function(fit, strata, seed = 1L, n_perm = 199L, alpha = 0.05, alpha_mem = 0.05, J = 6L,
                           resid_site = TRUE, min_cts = 2L, mem_rule = "holm", min_shared = 20L, rg = FALSE) {
  set.seed(seed)
  Z <- .prep_scores(fit, strata, resid_site)
  if (rg) Z <- lapply(Z, .rank_gauss)
  Z <- lapply(Z, .std_cols)
  out <- list(); id <- 0L
  for (j in seq_len(J)) {
    Z <- Z[vapply(Z, ncol, 1L) > 0]
    if (length(Z) < min_cts) break
    pairs <- .ct_pairs(Z, min_shared)
    if (!length(pairs)) break
    T <- .node_T(Z, pairs)
    lam <- .top_eig(T)$values[1]
    nl <- vapply(seq_len(n_perm), function(b) {
      Zp <- lapply(Z, function(S) .perm_rows(S, strata[rownames(S)]))
      .top_eig(.node_T(Zp, pairs))$values[1]
    }, 1)
    p <- (1 + sum(nl >= lam)) / (n_perm + 1)
    if (p >= alpha) break
    nct <- sub("\\|.*", "", rownames(T))
    tests <- do.call(rbind, lapply(names(Z), function(c) {
      keep <- nct != c
      u <- .top_eig(T[keep, keep, drop = FALSE])$vectors[, 1]; names(u) <- rownames(T)[keep]
      cbind(cell_type = c, .node_tests(Z[[c]], .loo_consensus(Z, u, c, rownames(Z[[c]]))))
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
