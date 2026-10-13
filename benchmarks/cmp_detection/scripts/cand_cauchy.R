## cand_cauchy.R -- Cauchy-combination node statistic (WS1 arm, winner-consistent prep).
## Replaces the winner's hard tau threshold with soft evidence accumulation. For node v (cell type c)
## and partner cell type c': m_{v,c'} = max_{u in c'} |T_vu| (same per-(node, partner-ct) max as
## cand_node_agg). Asymptotic p_{v,c'} = 2*pnorm(-m_{v,c'}) (Fisher-z is approximately N(0,1) under the
## null; this calibration is permutation-pool-independent, which is what keeps the subsequent permutation
## max-T exact). Node statistic:
##   S_v = mean_{c'} tan(pi*(0.5 - p_{v,c'}))          (Cauchy combination; p capped at [1e-15, 1-1e-15])
## Null: S_v recomputed under each of the n_perm within-strata permutations (.perm_Z), max over nodes per
##   permutation; fwer_p = (1 + #{maxS_perm >= S_v})/(n_perm+1); keep nodes with fwer_p < alpha.
## Grouping identical to the winner n_agg_noanc (anchors = FALSE): significant nodes linked by edges with
##   |T| >= link_thr, grouped with .cluster_edges over the all-nodes df.
## Prep identical to n_agg_noanc (.prep_rank_first(fit, strata, TRUE) + ambient drop at amb_thr = 0.3 via
##   ambient_scores(fit, fit$profiles)) so the comparison isolates the node statistic.
## Needs cand_edge.R (.prep_rank_first, .perm_Z, .edges_from), cand_nodeagg.R (.node_groups),
##   cand_spec.R (.ct_pairs), ambient.R (ambient_scores), methods.R (.perm_rows, .cluster_edges,
##   .empty_clusters).

## per-(node, partner-ct) max |T| -> asymptotic p -> Cauchy weights -> mean per node (named by "ct|prog")
.node_cauchy <- function(stat, G) {
  dt <- data.table::data.table(g = G$g, s = c(stat, stat))[, list(m = max(s)), by = "g"]
  p <- pmin(pmax(2 * stats::pnorm(-dt$m), 1e-15), 1 - 1e-15)
  dt[, w := tan(pi * (0.5 - p))]
  dt[, node := G$node[g]]
  out <- dt[, list(S = mean(w)), by = "node"]
  stats::setNames(out$S, as.character(out$node))
}

cand_cauchy <- function(fit, strata, seed = 1L, link_thr = 3, n_perm = 199L, alpha = 0.05,
                        min_shared = 20L, amb_thr = 0.3) {
  data.table::setDTthreads(1L)   # data.table's multithreaded grouping is pathologically slow on these
                                 # small per-perm aggregations (measured 8x slower than single-thread);
                                 # 1 thread matches the OMP_NUM_THREADS=1 intent of common.R. No effect
                                 # on results (grouping is deterministic).
  set.seed(seed)
  Z <- .prep_rank_first(fit, strata, TRUE)
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {   # same ambient drop as n_agg_noanc
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  pairs <- .ct_pairs(Z, min_shared)
  if (!length(pairs)) return(.empty_clusters())
  edges <- .edges_from(Z, pairs, FALSE)
  G <- .node_groups(edges)
  S <- .node_cauchy(edges$stat, G)
  mxS <- numeric(n_perm)
  for (b in seq_len(n_perm)) mxS[b] <- max(.node_cauchy(.edges_from(.perm_Z(Z, strata), pairs, FALSE)$stat, G))
  pS <- vapply(S, function(x) (1 + sum(mxS >= x)) / (n_perm + 1), 1)
  sig <- names(pS)[pS < alpha]
  keep <- G$ka %in% sig & G$kb %in% sig & edges$stat >= link_thr
  if (!any(keep)) return(.empty_clusters())
  E <- edges[keep, , drop = FALSE]
  E$fwer_p <- pmax(pS[G$ka[keep]], pS[G$kb[keep]])
  nodes <- do.call(rbind, lapply(names(Z), function(c) data.frame(cell_type = c, program_id = colnames(Z[[c]]), stringsAsFactors = FALSE)))
  .cluster_edges(E, nodes)
}
