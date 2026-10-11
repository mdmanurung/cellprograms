**Part 1 of 4 (covers your sections 1 and 3).** Nothing was stopped or changed for this. All jobs are still running. Everything below comes from read-only views of the benchmark scripts and result tables. One state note: my last action, appending FDR and NA-context tests to `test-coordination.R`, was cancelled by your pause. I checked afterwards and the tests were not written, so that file is unchanged.

**Parts of this knowledge transfer:**
1. **Part 1 (this reply):** section 1 (CMP callers: inputs, statistics, nulls, multiplicity control, clustering, membership, code) and section 3 (scoring: F1, node matching, upstream recall, the four false-call checks, code).
2. **Part 2:** section 2 (simulator: semi and param bases, how CMPs are planted, null and nullx).
3. **Part 3:** section 4 (pre-registration v1.1 and v1.2: decision rules, seed splits, deviations, adoption logic) and section 5 (screen results with failure modes).
4. **Part 4:** section 6 (real data: ambient, Ig/TCR, severity, site replication, covariates, K cap) and section 7 (operations).

---

# Section 1. CMP callers

## 1.0 Shared conventions

**Interface.** Every caller is `function(fit, strata, seed) -> clusters data.frame(cluster, cell_type, program_id, stat, fwer_p)`.

**Inputs.**
- `fit$scores[[ct]]`: donors × programs for each cell type (ct).
- `fit$scores_raw`: present only in within mode (unresidualized scores).
- `fit$loadings[[ct]]`: genes × programs.
- `fit$profiles`: mean CPM profile per cell type, used only by the ambient filter.
- `fit$ctx_var`: used only by C2.

**Strata.**
- Marginal mode: `strata` = site.
- Within mode (context scenario only): `strata` = site × context level. Context is binned to quintiles if it has more than 5 levels.

**Permutation unit.** Donor rows of each cell type's score matrix are permuted independently within strata. This breaks cross-cell-type coupling but keeps within-cell-type program correlation and site structure.

**Permutations and alpha.** n_perm = 199 in the benchmark (999 on real data). alpha = 0.05.

**Edge statistic.** For two nodes (ct_a, k) and (ct_b, l) on their shared donors, T = atanh(r) · √(n_shared − 3). Edges are computed only for cell-type pairs with at least 20 shared donors (`min_shared = 20`). The absolute value is used everywhere; the sign is stored but not used for calling.

**Grouping.** All callers except pa_cc use the repo's `.cluster_edges`: greedy single linkage over called edges, strongest first. A merge is skipped if it would put two nodes of the same cell type into one cluster. Clusters need at least 2 nodes.

**Helpers not pasted.** `.perm_within` and `sharing_spectrum` (used by pa_cc) are in your repo's `R/`. `.perm_rows`, `.cluster_edges`, `.empty_clusters`, `method_node_cor` and `method_pa_cc` are in your repo's `benchmarks/coordination/methods.R`. I used that file unchanged.

## 1.1 Summary table

| Caller (label) | Input | Statistic | Null | Multiplicity | Membership / grouping |
|---|---|---|---|---|---|
| **b_node_cor** (repo baseline) | raw scores (no site residualization, no rank transform) | edge T (Pearson) | within-site row permutation per ct | single-step max-T over all cross-ct edges | edges with FWER p < 0.05 → `.cluster_edges` |
| **b_pa_cc** (repo baseline) | scores, via `sharing_spectrum(space="scores")` | repo pairwise cell-type sharing test | repo's within-strata permutation (not restated here) | repo internal (not restated here) | connected components over *cell types*; `program_id = NA`, so cell-type-level clusters |
| **s_shuffled** (sanity) | scores, permuted once before testing | as `cand_edge` with Spearman-type ranks | max-T | max-T | should find nothing |
| **e_rf_maxT_gr_amb3** | rank-first scores + loadings/profiles for ambient filter | edge T on rank-Gauss scores | within-strata row permutation | single-step max-T over edges | significant edges → greedy |
| **n_agg_noanc** (v1 finalist, `main` in all v2 tables) | rank-first scores; loadings + profiles only for the ambient filter | node statistic A_v (below) | within-strata row permutation (199) | single-step max over nodes of A_v | significant nodes linked by edges with \|T\| ≥ 3 between two significant nodes → greedy |
| **n_agg_t2** | same | same as finalist, plus edge "anchors" | same | node max-T and edge max-T, each at alpha separately (not jointly controlled) | edges kept if edge FWER p < 0.05 **or** the node rule |
| **n_agg_minp_lt35** | same | anchors + link threshold 3.5 | same | single-step minP over the union of node and edge families | as n_agg_t2 |
| **c1_adapt** (v2 C1) | same | S_k(v) = sum of the k largest (m − tau)₊², k ∈ {1, 2, 4, ∞} | same | per-k max-T p → node min_k p, adjusted by the permutation distribution of per-permutation minP | as finalist |
| **c1_adapt_1a** | same | k ∈ {1, ∞} | same | same | same |
| **c2_ctx** (v2 C2) | as finalist plus `fit$ctx_var` | finalist node test at alpha/2, plus program-vs-context Fisher-z | node test as above; for context, context permuted within strata | Bonferroni split; context test is max-\|T\| over all programs | node CMPs as finalist; extra "context CMPs" = significant programs grouped by sign, strongest per ct, at least 2 cts, not already called |
| **c12** (`c12_ctx_adapt`) | C2 with C1 ks {1, 2, 4, ∞} | | | | |
| **j_r0 / j_r001** (MOFA arms) | joint-fit scores (MOFA factors projected per view) | as finalist (`cand_adaptive`, ks = ∞) | same | same | j_r001 first drops joint nodes with view R² < 0.001 |

The upstream arms in `screen_scorecard.csv` (main, K20, K30, K40, pnormal, featall, hvg4000) all use the frozen **n_agg_noanc** caller. Only the fit changes.

## 1.2 Finalist n_agg_noanc in words

1. **Rank-first preparation.** Each raw score column is rank-transformed to normal scores, then residualized on strata (site, or site × context). The order matters: ranking after residualizing let single-donor spikes shift a stratum mean, so the residual ranks ended up encoding strata.
2. **Ambient filter.** Nodes with ambient score > 0.3 are removed before testing (definition in 1.3).
3. **Edges.** All cross-cell-type node pairs with at least 20 shared donors.
4. **Node statistic.** For node v in cell type c, and each partner cell type c′, m_{v,c′} = max over nodes u in c′ of |T_vu|. Then A_v = Σ_{c′≠c} (m_{v,c′} − tau)₊², with tau = 2. This pools weak evidence over several partner cell types, so a program coupled to 3–4 cell types can be found even when no single edge passes edge max-T.
5. **Null.** 199 within-strata row permutations of every cell type. Each permutation records max_v A_v. Node p = (1 + #{perm max ≥ A_v}) / 200. Single-step max-T controls FWER over nodes, so for a fully null dataset P(any call) ≤ alpha.
6. **Membership.** An edge is kept if both ends have node p < 0.05 and |T| ≥ 3 (`link_thr`). Kept edges go to greedy single linkage with one node per cell type. Cluster `fwer_p` = max of the two node p's.
7. **"noanc".** With anchors on (n_agg_t2), any edge that is significant by edge max-T also enters. That family is not jointly controlled with the node family, which is why n_agg_t2 had null FWER 0.14 on the test set.

## 1.3 Core code (benchmark versions, verbatim)

```r
## ---- score preparation (cand_edge.R) ----
.rank_gauss <- function(S) {
  R <- apply(S, 2, function(v) stats::qnorm((rank(v, ties.method = "average") - 0.5) / length(v)))
  if (is.null(dim(R))) R <- matrix(R, nrow(S)); dimnames(R) <- dimnames(S); R
}
.prep_rank_first <- function(fit, strata, resid_site = TRUE, min_donors = 20L) {
  S0 <- if (!is.null(fit$scores_raw)) fit$scores_raw else fit$scores
  S0 <- S0[vapply(S0, function(S) !is.null(S) && ncol(S) > 0 && nrow(S) >= min_donors, TRUE)]
  lapply(S0, function(S) {
    R0 <- .rank_gauss(S); f <- factor(strata[rownames(S)])
    X <- if (resid_site && nlevels(f) > 1) stats::model.matrix(~f) else matrix(1, nrow(S), 1)
    R <- R0 - X %*% qr.solve(X, R0); dimnames(R) <- dimnames(S); R
  })
}
.perm_Z <- function(Z, strata) lapply(Z, function(S) .perm_rows(S, strata[rownames(S)]))  # .perm_rows: repo methods.R

## ---- pairs and edges (cand_spec.R / cand_edge.R) ----
.ct_pairs <- function(Z, min_shared = 20L) {
  cts <- names(Z); pairs <- list()
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    ids <- intersect(rownames(Z[[cts[i]]]), rownames(Z[[cts[j]]]))
    if (length(ids) >= min_shared) pairs[[length(pairs) + 1L]] <- list(a = cts[i], b = cts[j], ids = ids)
  }
  pairs
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

## ---- ambient score (ambient.R): max over other cts s of |cor(loading_k, log2(prof_s+1) - log2(prof_ct+1))|
ambient_scores <- function(fit, profiles) {
  lp <- lapply(profiles, function(p) log2(p + 1))
  do.call(rbind, lapply(fit$cell_types, function(ct) {
    W <- fit$loadings[[ct]]; if (is.null(W) || !ncol(W) || is.null(lp[[ct]])) return(NULL)
    best <- rep(0, ncol(W)); src <- rep(NA_character_, ncol(W))
    for (s in setdiff(names(lp), ct)) {
      g <- intersect(rownames(W), intersect(names(lp[[s]]), names(lp[[ct]])))
      if (length(g) < 50) next
      spec <- lp[[s]][g] - lp[[ct]][g]
      r <- abs(suppressWarnings(stats::cor(W[g, , drop = FALSE], spec))); r[!is.finite(r)] <- 0
      up <- r > best; best[up] <- r[up]; src[up] <- s
    }
    data.frame(cell_type = ct, program_id = colnames(W), amb = best, source = src, stringsAsFactors = FALSE)
  }))
}

## ---- node aggregation (cand_nodeagg.R): finalist = anchors = FALSE ----
.node_groups <- function(edges) {
  ka <- paste(edges$ct_a, edges$prog_a, sep = "|"); kb <- paste(edges$ct_b, edges$prog_b, sep = "|")
  ga <- paste(ka, edges$ct_b, sep = "@"); gb <- paste(kb, edges$ct_a, sep = "@")
  g <- factor(c(ga, gb)); nodes_of_g <- sub("@.*", "", levels(g))
  list(g = as.integer(g), ng = nlevels(g), node = factor(nodes_of_g), ka = ka, kb = kb)
}
.node_stat <- function(stat, G, tau) {      # A_v = sum over partner cts of (max |T| - tau)_+^2
  s <- c(stat, stat)
  m <- rep(-Inf, G$ng)
  dt <- data.table::data.table(g = G$g, s = s)[, list(m = max(s)), by = "g"]
  m[dt$g] <- dt$m
  v <- pmax(0, m - tau)^2
  rowsum(v, G$node, reorder = TRUE)[, 1]
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
  if (minp && anchors) {   # single-step minP over node and edge families
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
## registry entries
## n_agg_noanc     = cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = 199, anchors = FALSE)
## n_agg_t2        = cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3, n_perm = 199)          # anchors on
## n_agg_minp_lt35 = cand_node_agg(fit, strata, seed, tau = 2, link_thr = 3.5, n_perm = 199, minp = TRUE)
```

Edge-family callers (`cand_edge`; e_rf_maxT_gr_amb3 = rank-Gauss, max-T, greedy, amb_thr 0.3, rank_first TRUE):

```r
cand_edge <- function(fit, strata, seed = 1L, stat = "spearman", control = "maxT", q = 0.05, cluster = "greedy",
                      n_perm = 199L, alpha = 0.05, resid_site = TRUE, min_shared = 20L, resolution = 1,
                      shuffle = FALSE, amb_thr = Inf, rank_first = FALSE) {
  set.seed(seed)
  Z <- if (rank_first) .prep_rank_first(fit, strata, resid_site) else .transform_scores(.prep_scores(fit, strata, resid_site), stat)
  ## [same ambient filter block as cand_node_agg]
  if (shuffle) Z <- .perm_Z(Z, strata)            # s_shuffled: destroys all cross-ct coupling
  pairs <- .ct_pairs(Z, min_shared); if (!length(pairs)) return(.empty_clusters())
  edges <- .edges_from(Z, pairs, FALSE)
  if (control == "maxT") {
    mx <- vapply(seq_len(n_perm), function(b) max(.edges_from(.perm_Z(Z, strata), pairs, FALSE)$stat), 1)
    edges$fwer_p <- vapply(edges$stat, function(s) (1 + sum(mx >= s)) / (n_perm + 1), 1)
    edges <- edges[edges$fwer_p < alpha, , drop = FALSE]
  } ## (also "sdmaxT" Westfall-Young step-down and "bh" pooled-permutation BH branches; not in the final screens)
  if (!nrow(edges)) return(.empty_clusters())
  nodes <- do.call(rbind, lapply(names(Z), function(c) data.frame(cell_type = c, program_id = colnames(Z[[c]]), stringsAsFactors = FALSE)))
  if (cluster == "greedy") .cluster_edges(edges, nodes) else .leiden_group(edges, resolution)
}
## .prep_scores (non-rank-first): S - X qr.solve(X, S) with X = model.matrix(~site); .transform_scores applies .rank_gauss.
```

## 1.4 v2 callers (cand_v2.R, verbatim)

```r
## C1 node statistic for several k: S_k(v) = sum of the k largest (m_{v,c'} - tau)_+^2
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
.prep_Z_v2 <- function(fit, strata, amb_thr) {   # rank-first + ambient filter (identical to finalist)
  Z <- .prep_rank_first(fit, strata, TRUE)
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  Z
}
## ks = Inf reproduces n_agg_noanc exactly (verified identical on 6 datasets)
cand_adaptive <- function(fit, strata, seed = 1L, ks = c(1, 2, 4, Inf), tau = 2, link_thr = 3, n_perm = 199L,
                          alpha = 0.05, min_shared = 20L, amb_thr = 0.3, Z = NULL) {
  set.seed(seed)
  if (is.null(Z)) Z <- .prep_Z_v2(fit, strata, amb_thr)
  pairs <- .ct_pairs(Z, min_shared); if (!length(pairs)) return(.empty_clusters())
  edges <- .edges_from(Z, pairs, FALSE)
  G <- .node_groups(edges)
  S <- .node_topk(edges$stat, G, tau, ks)
  M <- matrix(0, n_perm, length(ks))
  for (b in seq_len(n_perm)) M[b, ] <- apply(.node_topk(.edges_from(.perm_Z(Z, strata), pairs, FALSE)$stat, G, tau, ks), 2, max)
  P <- sapply(seq_along(ks), function(j) vapply(S[, j], function(x) (1 + sum(M[, j] >= x)) / (n_perm + 1), 1))
  if (is.null(dim(P))) P <- matrix(P, ncol = length(ks))
  pmin_obs <- apply(P, 1, min)
  if (length(ks) > 1) {      # minP adjustment across k, using the same permutations
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
## C2: program-vs-context association (site-residualized context), T = Fisher z
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
  mx <- vapply(seq_len(n_perm), function(b) {      # context permuted within strata, max |T| over all programs
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
## Joint (MOFA) arms: drop joint nodes whose full-fit view R2 < thr, then cand_adaptive(ks = Inf)
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
REG_V2 <- list(
  c1_noanc_ref  = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = Inf, n_perm = NP),
  c1_adapt      = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = c(1, 2, 4, Inf), n_perm = NP),
  c1_adapt_1a   = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = c(1, Inf), n_perm = NP),
  c2_ctx        = function(fit, strata, seed) cand_ctx(fit, strata, seed, n_perm = NP, adaptive_ks = Inf),
  c12_ctx_adapt = function(fit, strata, seed) cand_ctx(fit, strata, seed, n_perm = NP, adaptive_ks = c(1, 2, 4, Inf)),
  j_r0          = function(fit, strata, seed) cand_adaptive(fit, strata, seed, ks = Inf, n_perm = NP),
  j_r001        = function(fit, strata, seed) cand_adaptive(.jfilter(fit, 0.001), strata, seed, ks = Inf, n_perm = NP)
)  # NP = 199
```

**How C2 gets a context variable in non-context scenarios** (from `run_cand_v2.R`). This lets the false-call checks see whether C2 inflates calls:

```r
ctx_for <- function(o) {
  if (!is.null(o$context)) return(o$context)          # true context (context scenario)
  ids <- rownames(o$meta); if (is.null(ids)) ids <- o$fit$observation_ids
  set.seed(as.integer((as.numeric(o$truth$seed) * 7919 + 17) %% 2147483647))
  c0 <- as.numeric(stats::runif(length(ids)) < 0.4); if (stats::sd(c0) == 0) c0[1] <- 1 - c0[1]
  stats::setNames(as.numeric(scale(c0)), ids)          # independent dummy binary, prevalence 0.4
}
## marginal mode: fit$ctx_var <- ctx_for(o); within mode: no ctx_var (C2 falls back to the node test at full alpha)
```

**Callers tried in v1 but not in the v2 screens** (registry names only; code available if you want it pasted):
- j_loo_mcca and j_node_spec: leave-one-cell-type-out MCCA or spectral membership tests, with Holm adjustment.
- j_flash.
- e_* variants: Pearson, winsorized, BH, Leiden, no-site.
- n_agg_* tuning variants: min_shared 15/30, link 2.5/3.5, amb 0.2/0.5, tau 2.5, no ambient filter, minP.
- e_rf_ext and n_agg_ext: extension step.
- b_subspace_sum and b_gene_match: repo baselines.

---

# Section 3. Scoring

## 3.1 Node-to-truth matching (`node_assign`)

**Candidate signals for a node.** A node is (cell type, program). Its candidates are the planted signals that include its cell type as a member. Their kinds are targets, decoys (ambient, genes_only), and context (a target in marginal mode, a decoy in within mode). A pseudo-signal **SITE** is added.

**Correlation scores.**
- For each candidate signal: |cor(node score, true donor latent z_{signal,ct})| on that cell type's donors.
- For SITE: √R² of lm(score ~ site).

**Assignment rule.** The node is assigned to the best-scoring candidate if that |cor| ≥ **0.4** (`NODE_THR`). Otherwise it is unassigned (NA). The assignment uses the scores that the caller saw in that mode; in within mode, those are the context-residualized scores. The sign of the correlation is ignored.

**CMP labels.** Each called CMP takes the majority signal of its nodes, which must cover more than 50% of its nodes. Otherwise it is "spurious". Label types are target, decoy, site, or spurious.

## 3.2 F1 per target signal (`eval_calls`)

For a target signal s with member cell types `mem`:
- For each called CMP c:
  - recall = (number of *distinct cell types* in c with a node assigned to s) / |mem|;
  - precision = (nodes of c assigned to s) / (nodes in c);
  - F1_c = harmonic mean of the two.
- **F1(s) = max over called CMPs of F1_c.** It is 0 if nothing is called.

The best-matching CMP is chosen per signal, and no one-to-one matching is enforced. Two signals could in principle both match the same CMP.

**Cell-type-level callers (pa_cc).** These have `program_id = NA`:
- recall counts member cell types in the cluster that have any node assigned to s;
- precision = |cluster ∩ mem| / cluster size.

This is more lenient than node-level matching, so it favours pa_cc.

## 3.3 Aggregates

- **Primary endpoint `membership_F1_pooled`.** Mean of F1(s) over all target-signal rows in the positive scenarios, marginal mode. Positive scenarios: cmp_distinct, cmp_mixed, context, ambient, lr, weak.
  - This is a mean over *signal rows*. Datasets with more targets weigh more.
  - In the tables, "F1" is this pooled number, and param_F1 / semi_F1 are the same per base.
- **Paired difference vs reference** (`v2_score.R`).
  - For each dataset (base, scenario, raw seed), take the mean F1 over its target signals.
  - Use only datasets common to all compared entries.
  - diff = candidate − reference per dataset, with a percentile CI from 2000 bootstrap resamples of datasets (`set.seed(1)`).
- **Breadth strata.** These are parsed from signal names: names containing "R20" → R20 (CMP spans ~20% of cell types, i.e. about 2), "R50" → R50, "CTX" → context, everything else → R100. `miss_<breadth>` = share of those signals with F1 exactly 0.
- **"Share of true member cell types with a matching program"** is `oracle` = upstream recall. For each target signal it is (number of distinct member cell types having ≥1 node assigned to s, anywhere in the fit) / |mem|. It does not depend on the caller and is a ceiling for recall. In the tables it appears as `upstream_recall` / `param_recall` / `semi_recall`.

## 3.4 The four false-call checks (dataset-level, binary, then averaged)

A dataset counts as a false call if it has *any* CMP of the relevant label type:

| Check | Datasets | Event counted |
|---|---|---|
| **Null FWER** (`*_null_fwer`) | null scenario: 30 dev seeds + 100 nullx seeds = 130 per base; in held-out, 50 | any CMP with label ≠ target. Null has no targets, so any CMP counts (spurious or site). Upper bound reported as a Clopper–Pearson 95% upper limit (`*_null_cp_hi`) |
| **Ambient decoy** | ambient scenario, marginal | any CMP labelled decoy (the leaked-ambient signal) |
| **Within-context false** | context scenario, within mode | any non-target CMP. In within mode the context signal is a decoy, so recovering it counts as false |
| **genes_only false** | genes_only scenario (shared gene loadings, independent donor scores) | any non-target CMP |

**Bounds.** The guardrail is ≤ 0.10 on each check, per base. `max_guardrail` = max over the four checks and both bases.

**Other diagnostics, not guardrails.** `site_any` (any CMP labelled site, any scenario) and `spurious_pos` (any spurious CMP in positive scenarios).

**Within-mode transforms used for evaluation:**

```r
within_strata <- function(strata, context) {
  if (is.null(context)) return(strata)
  cx <- context[names(strata)]
  lev <- if (length(unique(cx)) > 5) cut(cx, unique(stats::quantile(cx, 0:5 / 5)), include.lowest = TRUE) else cx
  stats::setNames(paste(strata, lev, sep = ":"), names(strata))
}
within_context <- function(fit, context) {   # residualize each ct's scores on context (+ intercept)
  if (is.null(context)) return(fit)
  fit$scores <- lapply(fit$scores, function(S) {
    if (!ncol(S)) return(S)
    x <- context[rownames(S)]; X <- cbind(1, x)
    R <- S - X %*% solve(crossprod(X), crossprod(X, S)); dimnames(R) <- dimnames(S); R
  })
  fit
}
## in within mode the caller also gets fit$scores_raw <- original scores, so rank-first re-residualizes on site x context
```

## 3.5 Scoring code (eval_cmp.R, verbatim core)

```r
NODE_THR <- 0.4
POS_SCEN <- c("cmp_distinct", "cmp_mixed", "context", "ambient", "lr", "weak")

.signal_table <- function(truth, strata, mode) {
  sg <- truth$signals
  typ <- vapply(sg, function(s) if (s$type == "context") (if (mode == "within") "decoy" else "target") else s$type, "")
  list(sg = sg, type = setNames(typ, names(sg)))
}
node_assign <- function(fit, truth, strata, mode = "marginal") {
  st <- .signal_table(truth, strata, mode); out <- list()
  for (ct in fit$cell_types) {
    S <- fit$scores[[ct]]; if (is.null(S) || !ncol(S)) next
    d <- rownames(S)
    cand <- names(st$sg)[vapply(st$sg, function(s) ct %in% s$members, TRUE)]
    C <- matrix(0, ncol(S), length(cand) + 1L, dimnames = list(colnames(S), c(cand, "SITE")))
    for (s in cand) { z <- st$sg[[s]]$z[[ct]][d]; C[, s] <- abs(suppressWarnings(stats::cor(S, z, use = "pairwise"))) }
    f <- factor(strata[d])
    C[, "SITE"] <- if (nlevels(f) > 1) sqrt(pmax(0, apply(S, 2, function(v) summary(stats::lm(v ~ f))$r.squared))) else 0
    C[is.na(C)] <- 0
    best <- colnames(C)[max.col(C, ties.method = "first")]; bv <- apply(C, 1, max)
    out[[ct]] <- data.frame(cell_type = ct, program_id = colnames(S), signal = ifelse(bv >= NODE_THR, best, NA),
                            cor = bv, stringsAsFactors = FALSE)
  }
  na <- do.call(rbind, out)
  na$type <- ifelse(is.na(na$signal), NA, ifelse(na$signal == "SITE", "site", st$type[na$signal]))
  na
}
eval_calls <- function(cl, fit, truth, strata, mode = "marginal") {
  na <- node_assign(fit, truth, strata, mode); st <- .signal_table(truth, strata, mode)
  tgt <- names(st$type)[st$type == "target"]
  ct_level <- nrow(cl) > 0 && all(is.na(cl$program_id))
  key <- paste(na$cell_type, na$program_id)
  cmps <- if (nrow(cl)) split(cl, cl$cluster) else list()
  lab <- vapply(cmps, function(c) {
    if (ct_level) {
      sigs <- unique(na$signal[!is.na(na$signal) & na$cell_type %in% c$cell_type])
      if (!length(sigs)) return("spurious")
      cov <- vapply(sigs, function(s) length(intersect(c$cell_type, na$cell_type[na$signal %in% s])), 1)
      return(if (max(cov) / nrow(c) > 0.5) sigs[which.max(cov)] else "spurious")
    }
    s <- na$signal[match(paste(c$cell_type, c$program_id), key)]
    tb <- table(s, useNA = "no")
    if (!length(tb) || max(tb) / nrow(c) <= 0.5) "spurious" else names(tb)[which.max(tb)]
  }, "")
  lab_type <- ifelse(lab == "spurious", "spurious", ifelse(lab == "SITE", "site", st$type[lab]))
  f1 <- vapply(tgt, function(s) {
    mem <- st$sg[[s]]$members
    if (!length(cmps)) return(0)
    max(vapply(cmps, function(c) {
      if (ct_level) {
        hit <- intersect(c$cell_type, intersect(mem, na$cell_type[na$signal %in% s]))
        rec <- length(hit) / length(mem); pre <- length(intersect(c$cell_type, mem)) / nrow(c)
      } else {
        s_n <- na$signal[match(paste(c$cell_type, c$program_id), key)]
        ok <- !is.na(s_n) & s_n == s
        rec <- length(unique(c$cell_type[ok])) / length(mem); pre <- mean(ok)
      }
      if (rec + pre == 0) 0 else 2 * rec * pre / (rec + pre)
    }, 1))
  }, 1)
  oracle <- vapply(tgt, function(s) length(unique(na$cell_type[na$signal %in% s])) / length(st$sg[[s]]$members), 1)
  list(f1 = f1, oracle = oracle, labels = lab, label_type = lab_type,
       summary = data.frame(n_cmp = length(cmps), n_false = sum(lab_type != "target"),
                            any_false = any(lab_type != "target"), any_decoy = any(lab_type == "decoy"),
                            any_site = any(lab_type == "site"), any_spurious = any(lab_type == "spurious"),
                            mean_f1 = if (length(f1)) mean(f1) else NA_real_,
                            mean_oracle = if (length(oracle)) mean(oracle) else NA_real_,
                            n_nodes_called = nrow(cl), stringsAsFactors = FALSE))
}
```

**Paired comparison and guardrails** (`v2_score.R`, core lines):

```r
seed_raw <- function(s) s %% 100000      # truth$seed carries a base/scenario offset
keyf <- function(g) unique(g[scenario %in% POS_SCEN & mode == "marginal", paste(base, scenario, sr)])
common <- Reduce(intersect, lapply(D, function(d) keyf(d$g)))      # identical datasets for all entries
brd <- function(sig) fifelse(grepl("R20", sig), "R20", fifelse(grepl("R50", sig), "R50", fifelse(sig == "CTX", "CTX", "R100")))
## per base b:  F1 = mean(f1), recall = mean(oracle), miss_<br> = mean(f1 == 0)
##   null: any_false over null-scenario dev rows + nullx rows; CP upper via binom.test
##   ambient_decoy = mean(any_decoy | ambient); within_false = mean(any_false | context, within); genes_only_false = mean(any_false | genes_only)
agg <- function(gg) gg[scenario %in% POS_SCEN & mode == "marginal" & paste(base, scenario, sr) %in% common, .(f1 = mean(f1)), by = .(base, scenario, sr)]
m <- merge(agg(ref$g), agg(d$g), by = c("base", "scenario", "sr"))
dd <- m$f1.y - m$f1.x; set.seed(1); bs <- replicate(2000, mean(sample(dd, replace = TRUE)))
diff_vs_ref <- mean(dd); diff_lo <- quantile(bs, .025); diff_hi <- quantile(bs, .975)
max_guardrail <- max(<all *_null_fwer, *_ambient_decoy, *_within_false, *_genes_only_false>)
```

## 3.6 Scoring caveats worth carrying into the repo

- **The 0.4 cut-off drives both numbers.** It defines both "matching program" and F1. A program that mixes two planted signals at |cor| 0.35 each counts as unassigned and lowers precision.
- **Ambient decoys count only when they are coherent.** Calling an ambient decoy is a false call only if it forms a CMP whose nodes are mostly that decoy.
- **The null check catches site artefacts.** Site artefacts in null data show up as "site" or "spurious" CMPs and count against FWER.
- **Gene loadings are not scored.** F1 checks only donor-score membership. A CMP whose member programs match the right donors but have wrong gene content would still score well.
- **No one-to-one matching.** Merging two true CMPs into one big CMP costs precision for both, but is not otherwise penalized.

Part 2 (simulator) comes next on request.
