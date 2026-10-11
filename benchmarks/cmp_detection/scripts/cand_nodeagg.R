## cand_nodeagg.R -- F3 variant: node-level evidence aggregation across partner cell types.
## Edge max-T is power-limited when cell types share few donors: a true program coupled to 3-4 cell types can have
## every single edge just below the max-T threshold. Node statistic pools evidence over partner cell types:
##   A_v = sum_{c' != ct(v)} max(0, max_{u in c'} |T_vu| - tau)^2
## FWER over nodes by permutation max (same within-strata row permutations as edges). Significant nodes are linked by
## edges with |T| >= link_thr (and by any max-T significant edge), then grouped greedily with one node per cell type.
## Under the global null P(any node significant) <= alpha, so dataset-level FWER is controlled for null datasets.
.node_groups <- function(edges) {
  ka <- paste(edges$ct_a, edges$prog_a, sep = "|"); kb <- paste(edges$ct_b, edges$prog_b, sep = "|")
  ga <- paste(ka, edges$ct_b, sep = "@"); gb <- paste(kb, edges$ct_a, sep = "@")
  g <- factor(c(ga, gb)); nodes_of_g <- sub("@.*", "", levels(g))
  list(g = as.integer(g), ng = nlevels(g), node = factor(nodes_of_g), ka = ka, kb = kb)
}
.node_stat <- function(stat, G, tau) {
  s <- c(stat, stat)
  m <- rep(-Inf, G$ng); m <- pmax(m, -Inf)
  dt <- data.table::data.table(g = G$g, s = s)[, list(m = max(s)), by = "g"]
  m[dt$g] <- dt$m
  v <- pmax(0, m - tau)^2
  out <- rowsum(v, G$node, reorder = TRUE)[, 1]
  out
}

cand_node_agg <- function(fit, strata, seed = 1L, tau = 2, link_thr = 3, n_perm = 199L, alpha = 0.05,
                          resid_site = TRUE, min_shared = 20L, amb_thr = 0.3, rank_first = TRUE, anchors = TRUE, minp = FALSE) {
  set.seed(seed)
  Z <- if (rank_first) .prep_rank_first(fit, strata, resid_site) else .transform_scores(.prep_scores(fit, strata, resid_site), "rankgauss")
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  pairs <- .ct_pairs(Z, min_shared)
  if (!length(pairs)) return(.empty_clusters())
  edges <- .edges_from(Z, pairs, FALSE)
  G <- .node_groups(edges)
  A <- .node_stat(edges$stat, G, tau)
  mxA <- numeric(n_perm); mxE <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    sb <- .edges_from(.perm_Z(Z, strata), pairs, FALSE)$stat
    mxE[b] <- max(sb); mxA[b] <- max(.node_stat(sb, G, tau))
  }
  pA <- vapply(A, function(x) (1 + sum(mxA >= x)) / (n_perm + 1), 1)
  edges$fwer_p <- vapply(edges$stat, function(s) (1 + sum(mxE >= s)) / (n_perm + 1), 1)
  if (minp && anchors) {   # union of node and edge families: single-step minP over both maxima (joint FWER <= alpha)
    pE_b <- vapply(mxE, function(s) sum(mxE >= s) / n_perm, 1); pA_b <- vapply(mxA, function(s) sum(mxA >= s) / n_perm, 1)
    mp <- pmin(pE_b, pA_b)
    adj <- function(p) vapply(p, function(x) (1 + sum(mp <= x)) / (n_perm + 1), 1)
    pA <- stats::setNames(adj(pA), names(pA)); edges$fwer_p <- adj(edges$fwer_p)
  }
  sig <- names(pA)[pA < alpha]
  anc <- if (anchors) edges$fwer_p < alpha else rep(FALSE, nrow(edges))
  keep <- anc | (G$ka %in% sig & G$kb %in% sig & edges$stat >= link_thr)
  if (!any(keep)) return(.empty_clusters())
  E <- edges[keep, , drop = FALSE]
  E$fwer_p <- pmin(E$fwer_p, pmax(pA[paste(E$ct_a, E$prog_a, sep = "|")], pA[paste(E$ct_b, E$prog_b, sep = "|")]), na.rm = TRUE)
  nodes <- do.call(rbind, lapply(names(Z), function(c) data.frame(cell_type = c, program_id = colnames(Z[[c]]), stringsAsFactors = FALSE)))
  .cluster_edges(E, nodes)
}
