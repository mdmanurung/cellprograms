## cand_edge.R -- F2/F3/F4 family: node-pair edge statistic -> error control -> grouping.
## stat:    pearson | spearman | rankgauss (Pearson on rank-based inverse-normal scores) | winsor (clip at 3 MAD)
## control: maxT (FWER over all cross-cell-type node pairs) | bh (pooled-permutation p-values, BH at q)
## cluster: greedy (repo single linkage w/ one-node-per-cell-type) | leiden (igraph, then best node per ct)
## cmp_null: if TRUE, a called CMP must also beat the permutation null of its own size-1 edge strength
.rank_gauss <- function(S) {
  R <- apply(S, 2, function(v) stats::qnorm((rank(v, ties.method = "average") - 0.5) / length(v)))
  if (is.null(dim(R))) R <- matrix(R, nrow(S)); dimnames(R) <- dimnames(S); R
}
.winsor <- function(S, k = 3) apply(S, 2, function(v) { m <- stats::median(v); s <- stats::mad(v) + 1e-12; pmin(pmax(v, m - k * s), m + k * s) })
.transform_scores <- function(Z, stat) {
  switch(stat, pearson = Z, spearman = lapply(Z, .rank_gauss), rankgauss = lapply(Z, .rank_gauss),
         winsor = lapply(Z, function(S) { W <- .winsor(S); dimnames(W) <- dimnames(S); W }), stop("stat"))
}
.pair_T <- function(A, B, ids, spearman = FALSE) {
  r <- suppressWarnings(stats::cor(A[ids, , drop = FALSE], B[ids, , drop = FALSE], method = if (spearman) "spearman" else "pearson"))
  r[!is.finite(r)] <- 0
  atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(length(ids) - 3)
}
.edges_from <- function(Z, pairs, spearman = FALSE) {
  do.call(rbind, lapply(pairs, function(p) {
    t <- .pair_T(Z[[p$a]], Z[[p$b]], p$ids, spearman)
    data.frame(ct_a = p$a, ct_b = p$b, prog_a = rownames(t)[row(t)], prog_b = colnames(t)[col(t)],
               stat = as.vector(abs(t)), sgn = as.vector(sign(t)), stringsAsFactors = FALSE)
  }))
}
## rank-first preparation: rank-gauss each raw score column, then residualize on strata (site, or site x context).
## Avoids the artefact where a single-donor spike shifts a stratum mean and the residual ranks then encode strata.
.prep_rank_first <- function(fit, strata, resid_site = TRUE, min_donors = 20L) {
  S0 <- if (!is.null(fit$scores_raw)) fit$scores_raw else fit$scores
  S0 <- S0[vapply(S0, function(S) !is.null(S) && ncol(S) > 0 && nrow(S) >= min_donors, TRUE)]
  lapply(S0, function(S) {
    R0 <- .rank_gauss(S); f <- factor(strata[rownames(S)])
    X <- if (resid_site && nlevels(f) > 1) stats::model.matrix(~f) else matrix(1, nrow(S), 1)
    R <- R0 - X %*% qr.solve(X, R0); dimnames(R) <- dimnames(S); R
  })
}
## rows of each score matrix get within-site permutations. For rank stats the ranks are permuted too (exact).
.perm_Z <- function(Z, strata) lapply(Z, function(S) .perm_rows(S, strata[rownames(S)]))

.leiden_group <- function(edges, resolution = 1) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph needed")
  ka <- paste(edges$ct_a, edges$prog_a, sep = "|"); kb <- paste(edges$ct_b, edges$prog_b, sep = "|")
  g <- igraph::graph_from_data_frame(data.frame(ka, kb, weight = edges$stat), directed = FALSE)
  set.seed(1)
  cm <- igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution_parameter = resolution)
  mem <- igraph::membership(cm); st <- igraph::strength(g)
  out <- list(); id <- 0L
  for (k in unique(mem)) {
    v <- names(mem)[mem == k]; ct <- sub("\\|.*", "", v)
    v <- v[order(-st[v])]; v <- v[!duplicated(sub("\\|.*", "", v))]   # one node per cell type, strongest
    if (length(v) < 2L) next
    id <- id + 1L
    e <- edges[ka %in% v & kb %in% v, , drop = FALSE]
    out[[id]] <- data.frame(cluster = id, cell_type = sub("\\|.*", "", v), program_id = sub("^[^|]*\\|", "", v),
                            stat = if (nrow(e)) max(e$stat) else NA, fwer_p = if (nrow(e)) min(e$fwer_p) else NA,
                            stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else .empty_clusters()
}

cand_edge <- function(fit, strata, seed = 1L, stat = "spearman", control = "maxT", q = 0.05, cluster = "greedy",
                      n_perm = 199L, alpha = 0.05, resid_site = TRUE, min_shared = 20L, resolution = 1,
                      shuffle = FALSE, amb_thr = Inf, rank_first = FALSE) {
  set.seed(seed)
  Z <- if (rank_first) .prep_rank_first(fit, strata, resid_site) else .transform_scores(.prep_scores(fit, strata, resid_site), stat)
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {   # F8: drop ambient-like nodes before testing
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  if (shuffle) Z <- .perm_Z(Z, strata)            # sanity candidate: destroys all cross-ct coupling
  pairs <- .ct_pairs(Z, min_shared)
  if (!length(pairs)) return(.empty_clusters())
  sp <- FALSE                                         # rankgauss/spearman handled by transform (Pearson on normal scores)
  edges <- .edges_from(Z, pairs, sp)
  if (control == "maxT") {
    mx <- vapply(seq_len(n_perm), function(b) max(.edges_from(.perm_Z(Z, strata), pairs, sp)$stat), 1)
    edges$fwer_p <- vapply(edges$stat, function(s) (1 + sum(mx >= s)) / (n_perm + 1), 1)
    edges <- edges[edges$fwer_p < alpha, , drop = FALSE]
  } else if (control == "sdmaxT") {   # Westfall-Young step-down max-T (strong FWER control, more power than single-step)
    P <- vapply(seq_len(n_perm), function(b) .edges_from(.perm_Z(Z, strata), pairs, sp)$stat, numeric(nrow(edges)))
    o <- order(edges$stat, decreasing = TRUE); Po <- P[o, , drop = FALSE]
    ## successive maxima from the bottom: M[j, b] = max_{k >= j} Po[k, b]
    M <- apply(Po, 2, function(v) rev(cummax(rev(v))))
    if (is.null(dim(M))) M <- matrix(M, nrow = 1)
    pj <- vapply(seq_along(o), function(j) (1 + sum(M[j, ] >= edges$stat[o[j]])) / (n_perm + 1), 1)
    edges$fwer_p[o] <- cummax(pj)
    edges <- edges[edges$fwer_p < alpha, , drop = FALSE]
  } else {
    nb <- max(20L, n_perm %/% 5L)
    nullv <- unlist(lapply(seq_len(nb), function(b) .edges_from(.perm_Z(Z, strata), pairs, sp)$stat))
    nullv <- sort(nullv)
    edges$fwer_p <- (1 + length(nullv) - findInterval(edges$stat, nullv, left.open = TRUE)) / (length(nullv) + 1)
    edges$fwer_p <- pmin(1, edges$fwer_p)
    edges <- edges[stats::p.adjust(edges$fwer_p, "BH") < q, , drop = FALSE]
  }
  if (!nrow(edges)) return(.empty_clusters())
  nodes <- do.call(rbind, lapply(names(Z), function(c) data.frame(cell_type = c, program_id = colnames(Z[[c]]), stringsAsFactors = FALSE)))
  if (cluster == "greedy") .cluster_edges(edges, nodes) else .leiden_group(edges, resolution)
}
