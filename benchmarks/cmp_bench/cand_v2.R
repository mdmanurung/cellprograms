## cand_v2.R -- v2 caller arms (PLAN v2, option 3). Requires registry.R sourced first.
## C1 breadth-adaptive node test: per node v, m_{v,c'} = max edge |T| to partner cell type c';
##   S_k(v) = sum of the k largest (m - tau)_+^2, k in ks (Inf = all partners; k = 1 targets 2-cell-type CMPs).
##   Per-k single-step max-T p over nodes; node p = min_k p_k, adjusted by the permutation distribution of the
##   per-permutation minP (same 199 within-strata permutations) -> FWER over nodes and k jointly <= alpha.
## C2 context-supervised add-on: node test at alpha/2 + context-association test at alpha/2 (Bonferroni split).
##   Context test: signed Fisher-z of (site-residualized rank-gauss score, site-residualized context) per program,
##   context permuted within strata, max-|T| over programs. Significant programs grouped by sign (one per cell type,
##   strongest); a group with >= 2 cell types and not already in a node-test CMP is a "context-driven" CMP.

.node_topk <- function(stat, G, tau, ks) {
  dt <- data.table::data.table(g = G$g, s = c(stat, stat))[, list(m = max(s)), by = "g"]
  dt[, node := G$node[g]][, v := pmax(0, m - tau)^2]
  data.table::setorder(dt, node, -v)
  dt[, r := seq_len(.N), by = "node"][, cs := cumsum(v), by = "node"]
  nodes <- levels(G$node)
  out <- vapply(ks, function(k) {
    kk <- if (is.finite(k)) k else .Machine$integer.max
    x <- dt[r <= kk, list(S = max(cs)), by = "node"]
    s <- stats::setNames(numeric(length(nodes)), nodes); s[as.character(x$node)] <- x$S; s
  }, numeric(length(nodes)))
  if (is.null(dim(out))) out <- matrix(out, ncol = length(ks))
  rownames(out) <- nodes; out
}

.prep_Z_v2 <- function(fit, strata, amb_thr) {
  Z <- .prep_rank_first(fit, strata, TRUE)
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  Z
}

## node test returning clusters; ks = Inf reproduces n_agg_noanc exactly (S_Inf = A_v)
cand_adaptive <- function(fit, strata, seed = 1L, ks = c(1, 2, 4, Inf), tau = 2, link_thr = 3, n_perm = 199L,
                          alpha = 0.05, min_shared = 20L, amb_thr = 0.3, Z = NULL) {
  set.seed(seed)
  if (is.null(Z)) Z <- .prep_Z_v2(fit, strata, amb_thr)
  pairs <- .ct_pairs(Z, min_shared)
  if (!length(pairs)) return(.empty_clusters())
  edges <- .edges_from(Z, pairs, FALSE)
  G <- .node_groups(edges)
  S <- .node_topk(edges$stat, G, tau, ks)
  M <- matrix(0, n_perm, length(ks))
  for (b in seq_len(n_perm)) M[b, ] <- apply(.node_topk(.edges_from(.perm_Z(Z, strata), pairs, FALSE)$stat, G, tau, ks), 2, max)
  P <- sapply(seq_along(ks), function(j) vapply(S[, j], function(x) (1 + sum(M[, j] >= x)) / (n_perm + 1), 1))
  if (is.null(dim(P))) P <- matrix(P, ncol = length(ks))
  pmin_obs <- apply(P, 1, min)
  if (length(ks) > 1) {
    Pb <- sapply(seq_along(ks), function(j) vapply(M[, j], function(x) sum(M[, j] >= x) / n_perm, 1))
    q <- apply(Pb, 1, min)
    pA <- vapply(pmin_obs, function(x) (1 + sum(q <= x)) / (n_perm + 1), 1)
  } else pA <- pmin_obs
  names(pA) <- rownames(S)
  sig <- names(pA)[pA < alpha]
  keep <- G$ka %in% sig & G$kb %in% sig & edges$stat >= link_thr
  if (!any(keep)) return(.empty_clusters())
  E <- edges[keep, , drop = FALSE]
  E$fwer_p <- pmax(pA[G$ka[keep]], pA[G$kb[keep]])
  nodes <- do.call(rbind, lapply(names(Z), function(c) data.frame(cell_type = c, program_id = colnames(Z[[c]]), stringsAsFactors = FALSE)))
  .cluster_edges(E, nodes)
}

.ctx_assoc_T <- function(Z, ctx, strata) {
  do.call(rbind, lapply(names(Z), function(c) {
    S <- Z[[c]]; ids <- rownames(S); x <- ctx[ids]; ok <- is.finite(x)
    if (sum(ok) < 10 || stats::sd(x[ok]) == 0) return(NULL)
    f <- factor(strata[ids][ok]); X <- if (nlevels(f) > 1) stats::model.matrix(~f) else matrix(1, sum(ok), 1)
    xr <- x[ok] - X %*% qr.solve(X, x[ok])
    r <- suppressWarnings(stats::cor(S[ok, , drop = FALSE], xr)); r[!is.finite(r)] <- 0
    data.frame(cell_type = c, program_id = colnames(S), T = atanh(pmin(pmax(as.vector(r), -0.999999), 0.999999)) * sqrt(sum(ok) - 3),
               stringsAsFactors = FALSE)
  }))
}

cand_ctx <- function(fit, strata, seed = 1L, alpha = 0.05, n_perm = 199L, amb_thr = 0.3, adaptive_ks = Inf, ...) {
  ctx <- fit$ctx_var
  if (is.null(ctx)) return(cand_adaptive(fit, strata, seed, ks = adaptive_ks, alpha = alpha, n_perm = n_perm, amb_thr = amb_thr, ...))
  Z <- .prep_Z_v2(fit, strata, amb_thr)
  cl <- cand_adaptive(fit, strata, seed, ks = adaptive_ks, alpha = alpha / 2, n_perm = n_perm, amb_thr = amb_thr, Z = Z, ...)
  set.seed(seed + 7L)
  obs <- .ctx_assoc_T(Z, ctx, strata)
  if (is.null(obs) || !nrow(obs)) return(cl)
  ids <- names(ctx); st <- strata[ids]
  mx <- vapply(seq_len(n_perm), function(b) {
    cp <- ctx; for (s in unique(st)) { w <- which(st == s); cp[w] <- ctx[w][sample.int(length(w))] }
    max(abs(.ctx_assoc_T(Z, stats::setNames(cp, ids), strata)$T))
  }, 1)
  obs$p <- vapply(abs(obs$T), function(x) (1 + sum(mx >= x)) / (n_perm + 1), 1)
  sig <- obs[obs$p < alpha / 2, , drop = FALSE]
  if (nrow(cl)) sig <- sig[!paste(sig$cell_type, sig$program_id) %in% paste(cl$cell_type, cl$program_id), , drop = FALSE]
  id0 <- if (nrow(cl)) max(cl$cluster) else 0L
  add <- list()
  for (sg in c(1, -1)) {
    g <- sig[sign(sig$T) == sg, , drop = FALSE]
    g <- g[order(-abs(g$T)), , drop = FALSE]; g <- g[!duplicated(g$cell_type), , drop = FALSE]
    if (nrow(g) >= 2) { id0 <- id0 + 1L
      add[[length(add) + 1]] <- data.frame(cluster = id0, cell_type = g$cell_type, program_id = g$program_id,
                                           stat = max(abs(g$T)), fwer_p = max(g$p), stringsAsFactors = FALSE) }
  }
  if (length(add)) rbind(cl, do.call(rbind, add)) else cl
}

REG_V2 <- list(
  c1_noanc_ref = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = Inf, n_perm = NP),
  c1_adapt     = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = c(1, 2, 4, Inf), n_perm = NP),
  c1_adapt_1a  = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = c(1, Inf), n_perm = NP),
  c2_ctx       = function(fit, strata, seed) cand_ctx(fit, strata, seed, n_perm = NP, adaptive_ks = Inf),
  c12_ctx_adapt = function(fit, strata, seed) cand_ctx(fit, strata, seed, n_perm = NP, adaptive_ks = c(1, 2, 4, Inf))
)

## joint-arm node filter: drop joint (factor, view) nodes whose full-fit view R2 < thr (per-ct programs untouched)
.jfilter <- function(fit, thr) {
  r2 <- fit$joint$node_r2; if (is.null(r2)) return(fit)
  for (ct in names(fit$scores)) {
    S <- fit$scores[[ct]]; if (is.null(S) || !ncol(S)) next
    drop <- colnames(S) %in% names(r2)[r2 < thr]
    fit$scores[[ct]] <- S[, !drop, drop = FALSE]
    if (!is.null(fit$scores_raw[[ct]])) fit$scores_raw[[ct]] <- fit$scores_raw[[ct]][, !drop, drop = FALSE]
    if (!is.null(fit$loadings[[ct]])) fit$loadings[[ct]] <- fit$loadings[[ct]][, setdiff(colnames(fit$loadings[[ct]]), colnames(S)[drop]), drop = FALSE]
  }
  fit
}
REG_V2$j_r001 <- function(fit, strata, seed) cand_adaptive(.jfilter(fit, 0.001), strata, seed, ks = Inf, n_perm = NP)
REG_V2$j_r005 <- function(fit, strata, seed) cand_adaptive(.jfilter(fit, 0.005), strata, seed, ks = Inf, n_perm = NP)
REG_V2$j_r0   <- function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = Inf, n_perm = NP)

## S: subspace_sum joint-subspace arm (reviewed prototype on winner prep; cand_subspace.R)
source("/mnt/shared-workspace/cmp/bench/cand_subspace.R")
REG_V2$s_subspace      <- function(fit, strata, seed) cand_subspace(fit, strata, seed, gate = "maxT", n_perm = NP)
REG_V2$s_subspace_ct95 <- function(fit, strata, seed) cand_subspace(fit, strata, seed, gate = "ct95", n_perm = NP)

## WS1: two-stage pair-gate -> node test (cand_twostage.R) and Cauchy-combination node statistic
## (cand_cauchy.R); both on the identical n_agg_noanc prep so the comparison isolates the node statistic.
source("/mnt/shared-workspace/cmp/bench/cand_twostage.R")
source("/mnt/shared-workspace/cmp/bench/cand_cauchy.R")
REG_V2$n_2stage <- function(fit, strata, seed) cand_twostage(fit, strata, seed, n_perm = NP)
REG_V2$n_cauchy  <- function(fit, strata, seed) cand_cauchy(fit, strata, seed, n_perm = NP)

## WS2.1: score-uncertainty arm (cand_nodeagg_wt.R) -- EB shrinkage of donor scores by flashier
## L_psd (canonicalized), then the exact n_agg_noanc pipeline; marginal mode only (within mode
## falls back to the unmodified winner). Isolates uncertainty propagation as a caller arm.
source("/mnt/shared-workspace/cmp/bench/cand_nodeagg_wt.R")
REG_V2$n_agg_eb <- function(fit, strata, seed) cand_nodeagg_eb(fit, strata, seed, n_perm = NP)
