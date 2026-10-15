## cand_twostage.R -- two-stage pair-gate -> node test (WS1 arm, winner-consistent prep).
## Stage 1 (pair gate): all cross-cell-type edge Fisher-z statistics as in cand_node_agg
##   (.edges_from over the shared-donor pairs from .ct_pairs); keep edges passing a LIBERAL raw gate
##   |T| >= gate_t (default 2.0).
## Stage 2 (node test): winner node statistic A_v = sum_{c' != ct(v)} max(0, max_{u in c'} |T_vu| - tau)^2
##   computed ONLY over surviving edges (node groups rebuilt from the gated edge set).
## Error control: the ENTIRE two-stage procedure (gate + node statistic) is recomputed inside each of the
##   n_perm within-strata permutations (.perm_Z), giving the max-A null distribution; per-node
##   fwer_p = (1 + #{maxA_perm >= A_v})/(n_perm+1); keep nodes with fwer_p < alpha. The gate is
##   data-dependent, so it must be inside the permutation loop for validity.
## Grouping: significant nodes linked by (surviving) edges with |T| >= link_thr, grouped with
##   .cluster_edges over the all-nodes df, exactly as the winner n_agg_noanc (anchors = FALSE).
## Prep identical to n_agg_noanc (.prep_rank_first(fit, strata, TRUE) + ambient drop at amb_thr = 0.3 via
##   ambient_scores(fit, fit$profiles)) so the comparison isolates the node statistic.
## Needs cand_edge.R (.prep_rank_first, .perm_Z, .edges_from), cand_nodeagg.R (.node_groups, .node_stat),
##   cand_spec.R (.ct_pairs), ambient.R (ambient_scores), methods.R (.perm_rows, .cluster_edges,
##   .empty_clusters).

cand_twostage <- function(fit, strata, seed = 1L, tau = 2, gate_t = 2, link_thr = 3, n_perm = 199L,
                          alpha = 0.05, min_shared = 20L, amb_thr = 0.3) {
  data.table::setDTthreads(1L)   # data.table's multithreaded grouping is pathologically slow on these
                                 # small per-perm aggregations (measured 16x slower than single-thread);
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
  surv <- edges$stat >= gate_t                              # stage 1: liberal pair gate on raw |T|
  if (!any(surv)) return(.empty_clusters())
  Es <- edges[surv, , drop = FALSE]
  G <- .node_groups(Es)
  A <- .node_stat(Es$stat, G, tau)                          # stage 2: node statistic over surviving edges
  mxA <- numeric(n_perm)
  for (b in seq_len(n_perm)) {                              # whole two-stage procedure inside the loop
    eb <- .edges_from(.perm_Z(Z, strata), pairs, FALSE)
    kb <- eb$stat >= gate_t
    mxA[b] <- if (any(kb)) max(.node_stat(eb$stat[kb], .node_groups(eb[kb, , drop = FALSE]), tau)) else 0
  }
  pA <- vapply(A, function(x) (1 + sum(mxA >= x)) / (n_perm + 1), 1)
  sig <- names(pA)[pA < alpha]
  keep <- G$ka %in% sig & G$kb %in% sig & Es$stat >= link_thr
  if (!any(keep)) return(.empty_clusters())
  E <- Es[keep, , drop = FALSE]
  E$fwer_p <- pmax(pA[G$ka[keep]], pA[G$kb[keep]])
  nodes <- do.call(rbind, lapply(names(Z), function(c) data.frame(cell_type = c, program_id = colnames(Z[[c]]), stringsAsFactors = FALSE)))
  .cluster_edges(E, nodes)
}
