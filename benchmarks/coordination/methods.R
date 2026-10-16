## methods.R -- candidate methods for recurring programs across cell types.
## Each takes a canonicalized fit (fit$scores, fit$loadings, fit$cell_types) and returns
##   list(method, clusters = data.frame(cluster, cell_type, program_id, stat, fwer_p), n_perm)
## program_id is NA for pa_cc (cell-type-level clusters). Needs R/ sourced (.perm_within, .subspace, sharing_spectrum).

## ---- shared helpers -------------------------------------------------------
.nodes <- function(fit) {
  do.call(rbind, lapply(fit$cell_types, function(ct) {
    K <- ncol(fit$scores[[ct]])
    if (K == 0L) return(NULL)
    data.frame(cell_type = ct, program_id = colnames(fit$scores[[ct]]), stringsAsFactors = FALSE)
  }))
}

## permute the rows of M within strata (rows are donors of one cell type). Row names keep the
## ORIGINAL donor order: callers that re-index by donor name (node_cor) must see the permuted
## values under the original names, else the permutation is undone and the null equals the data.
.perm_rows <- function(M, strata_ct) {
  P <- M[.perm_within(nrow(M), strata_ct), , drop = FALSE]
  rownames(P) <- rownames(M)
  P
}

## Greedy single-linkage on called edges, strongest first; a merge that would put two
## nodes of the same cell type in one cluster is skipped. min_density > 0 drops clusters
## whose called edges cover fewer than that fraction of the possible cross-cell-type
## member pairs (a single-linkage chain A-B-C with a weak A-C passes the edge test but
## fails the density filter).
.cluster_edges <- function(edges, nodes, min_density = 0) {
  key <- paste(nodes$cell_type, nodes$program_id, sep = "|")
  cl <- setNames(seq_along(key), key)
  members <- setNames(as.list(key), seq_along(key))
  edges <- edges[order(-edges$stat), , drop = FALSE]
  used <- list()
  for (i in seq_len(nrow(edges))) {
    ka <- paste(edges$ct_a[i], edges$prog_a[i], sep = "|"); kb <- paste(edges$ct_b[i], edges$prog_b[i], sep = "|")
    ca <- cl[[ka]]; cb <- cl[[kb]]
    if (ca == cb) { used[[as.character(ca)]] <- c(used[[as.character(ca)]], i); next }
    cta <- sub("\\|.*", "", members[[as.character(ca)]]); ctb <- sub("\\|.*", "", members[[as.character(cb)]])
    if (length(intersect(cta, ctb))) next
    members[[as.character(ca)]] <- c(members[[as.character(ca)]], members[[as.character(cb)]])
    cl[members[[as.character(cb)]]] <- ca
    members[[as.character(cb)]] <- NULL
    used[[as.character(ca)]] <- c(used[[as.character(ca)]], used[[as.character(cb)]], i)
    used[[as.character(cb)]] <- NULL
  }
  out <- list(); id <- 0L
  for (nm in names(members)) {
    m <- members[[nm]]
    if (length(m) < 2L) next
    e <- edges[used[[nm]], , drop = FALSE]
    if (min_density > 0) {
      ctm <- sub("\\|.*", "", m)
      tab <- table(ctm)
      n_pairs <- choose(length(m), 2L) - sum(tab * (tab - 1L) / 2)
      if (nrow(e) / max(n_pairs, 1L) < min_density) next
    }
    id <- id + 1L
    out[[id]] <- data.frame(cluster = id, cell_type = sub("\\|.*", "", m), program_id = sub("^[^|]*\\|", "", m),
                            stat = max(e$stat), fwer_p = min(e$fwer_p), stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else .empty_clusters()
}
.empty_clusters <- function() data.frame(cluster = integer(), cell_type = character(), program_id = character(),
                                         stat = numeric(), fwer_p = numeric(), stringsAsFactors = FALSE)

## all cross-cell-type node pairs: statistic function applied per cell-type pair, max-T null
.maxT_edges <- function(pairs, stat_fun, perm_fun, n_perm, alpha) {
  obs <- lapply(pairs, stat_fun)
  edges <- do.call(rbind, lapply(seq_along(pairs), function(i) {
    o <- obs[[i]]; if (is.null(o)) return(NULL)
    data.frame(ct_a = pairs[[i]]$a, ct_b = pairs[[i]]$b,
               prog_a = rownames(o)[row(o)], prog_b = colnames(o)[col(o)], stat = as.vector(abs(o)),
               stringsAsFactors = FALSE)
  }))
  if (is.null(edges) || !nrow(edges)) return(edges)
  mx <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    pm <- perm_fun()
    mx[b] <- max(vapply(seq_along(pairs), function(i) {
      if (is.null(obs[[i]])) return(0)
      max(abs(stat_fun(pairs[[i]], pm)))
    }, numeric(1L)))
  }
  edges$fwer_p <- vapply(edges$stat, function(s) (1 + sum(mx >= s)) / (n_perm + 1), numeric(1L))
  edges[edges$fwer_p < alpha, , drop = FALSE]
}

## ---- node_cor: studentized score correlation, max-T over node pairs --------
method_node_cor <- function(fit, strata, n_perm = 999L, seed = 1L, alpha = 0.05, min_shared = 20L,
                            min_density = 0) {
  set.seed(seed)
  nodes <- .nodes(fit)
  cts <- unique(nodes$cell_type)
  Z <- fit$scores[cts]
  pairs <- list()
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    a <- cts[i]; b <- cts[j]
    ids <- intersect(rownames(Z[[a]]), rownames(Z[[b]]))
    if (length(ids) >= min_shared) pairs[[length(pairs) + 1L]] <- list(a = a, b = b, ids = ids)
  }
  stat_fun <- function(p, pm = NULL) {
    Za <- if (is.null(pm)) Z[[p$a]] else pm[[p$a]]; Zb <- if (is.null(pm)) Z[[p$b]] else pm[[p$b]]
    r <- suppressWarnings(stats::cor(Za[p$ids, , drop = FALSE], Zb[p$ids, , drop = FALSE]))
    r[!is.finite(r)] <- 0
    atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(length(p$ids) - 3)
  }
  perm_fun <- function() lapply(Z, function(M) .perm_rows(M, strata[rownames(M)]))
  edges <- .maxT_edges(pairs, stat_fun, perm_fun, n_perm, alpha)
  list(method = "node_cor", n_perm = n_perm,
       clusters = if (is.null(edges) || !nrow(edges)) .empty_clusters() else .cluster_edges(edges, nodes, min_density))
}

## ---- gene_match: loading cosine over shared genes, gene-label permutation, max-T ----
method_gene_match <- function(fit, n_perm = 999L, seed = 1L, alpha = 0.05, min_genes = 100L) {
  set.seed(seed)
  nodes <- .nodes(fit)
  cts <- unique(nodes$cell_type)
  W <- fit$loadings[cts]
  pairs <- list()
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    a <- cts[i]; b <- cts[j]
    g <- intersect(rownames(W[[a]]), rownames(W[[b]]))
    if (length(g) >= min_genes) pairs[[length(pairs) + 1L]] <- list(a = a, b = b, g = g)
  }
  unit <- function(M) { n <- sqrt(colSums(M^2)); n[n == 0] <- 1; sweep(M, 2, n, "/") }
  stat_fun <- function(p, pm = NULL) {
    Wa <- unit(W[[p$a]][p$g, , drop = FALSE])
    Wb <- if (is.null(pm)) unit(W[[p$b]][p$g, , drop = FALSE]) else unit(W[[p$b]][pm[[paste(p$a, p$b)]], , drop = FALSE])
    crossprod(Wa, Wb) * sqrt(length(p$g))
  }
  perm_fun <- local({
    pp <- pairs
    function() setNames(lapply(pp, function(p) match(sample(p$g), rownames(W[[p$b]]))), vapply(pp, function(p) paste(p$a, p$b), ""))
  })
  edges <- .maxT_edges(pairs, stat_fun, perm_fun, n_perm, alpha)
  list(method = "gene_match", n_perm = n_perm,
       clusters = if (is.null(edges) || !nrow(edges)) .empty_clusters() else .cluster_edges(edges, nodes))
}

## ---- pa_cc: existing pairwise test, connected components over cell types ----
method_pa_cc <- function(fit, strata, n_perm = 999L, seed = 1L) {
  ss <- sharing_spectrum(fit, space = "scores", n_perm = n_perm, seed = seed, strata = strata)$summary
  called <- ss[ss$shared, , drop = FALSE]
  cts <- fit$cell_types[vapply(fit$scores, ncol, 1L) > 0L]
  comp <- setNames(seq_along(cts), cts)
  for (i in seq_len(nrow(called))) {
    p <- strsplit(called$pair[i], " vs ", fixed = TRUE)[[1]]
    a <- comp[[p[1]]]; b <- comp[[p[2]]]
    if (a != b) comp[comp == b] <- a
  }
  g <- split(names(comp), comp); g <- g[lengths(g) >= 2L]
  cl <- if (!length(g)) .empty_clusters() else do.call(rbind, lapply(seq_along(g), function(i)
    data.frame(cluster = i, cell_type = g[[i]], program_id = NA_character_, stat = NA_real_, fwer_p = NA_real_,
               stringsAsFactors = FALSE)))
  list(method = "pa_cc", n_perm = n_perm, clusters = cl)
}

## ---- subspace_sum: sum of cell-type score projections, generalized eigenproblem ----
## M = sum_c U_c U_c' (zero rows for donors a cell type lacks); M v = lambda D v with
## D = diag(#cell types observing each donor), so lambda ~ fraction of cell types holding the direction.
## Sequential test with deflation: direction j is tested against a null built from the scores with
## directions 1..j-1 projected out of every cell type, then rows permuted within strata. Stop at the
## first p >= alpha (fixed sequence, FWER <= alpha). Using the all-permuted null for every j (the old
## cummax + break) ignores the signal already found and is not a valid step-down.
method_subspace_sum <- function(fit, strata, n_perm = 999L, seed = 1L, alpha = 0.05, J = 15L, min_shared = 20L) {
  set.seed(seed)
  nodes <- .nodes(fit)
  cts <- unique(nodes$cell_type)
  ids <- sort(unique(unlist(lapply(fit$scores[cts], rownames))))
  nu <- length(ids)
  Ub <- lapply(cts, function(ct) .subspace(fit$scores[[ct]])$U); names(Ub) <- cts
  for (ct in cts) rownames(Ub[[ct]]) <- rownames(fit$scores[[ct]])
  obs <- lapply(cts, function(ct) match(rownames(fit$scores[[ct]]), ids)); names(obs) <- cts
  Dd <- numeric(nu); for (ct in cts) Dd[obs[[ct]]] <- Dd[obs[[ct]]] + 1
  Dh <- 1 / sqrt(Dd)
  embed <- function(U, ct) { E <- matrix(0, nu, ncol(U)); E[obs[[ct]], ] <- U; E }
  deflate <- function(Us, V) {                      # remove span(V) (donor space) from each cell type
    if (!ncol(V)) return(Us)
    Map(function(U, ct) { Q <- qr.Q(qr(V[obs[[ct]], , drop = FALSE])); U - Q %*% crossprod(Q, U) }, Us, cts)
  }
  top <- function(Us) {                             # top direction, its value and per-cell-type contribution
    M <- Reduce(`+`, lapply(Map(embed, Us, cts), tcrossprod))
    ev <- eigen(Dh * M * rep(Dh, each = nu), symmetric = TRUE)
    v <- ev$vectors[, 1L] * Dh                      # back to donor space
    contrib <- vapply(seq_along(cts), function(k) {
      o <- obs[[k]]; sum(crossprod(Us[[k]], v[o])^2) / max(sum(v[o]^2), 1e-12)
    }, numeric(1L))
    list(value = ev$values[1L], v = v, contrib = contrib)
  }
  out <- list(); id <- 0L
  V <- matrix(0, nu, 0L)
  for (j in seq_len(J)) {
    Ud <- deflate(Ub, V)
    o <- top(Ud)
    nl <- numeric(n_perm); nc <- matrix(NA_real_, n_perm, length(cts))
    for (b in seq_len(n_perm)) {
      s <- top(Map(function(U, ct) .perm_rows(U, strata[rownames(fit$scores[[ct]])]), Ud, cts))
      nl[b] <- s$value; nc[b, ] <- s$contrib
    }
    p <- (1 + sum(nl >= o$value)) / (n_perm + 1)
    if (p >= alpha) break
    V <- cbind(V, o$v)
    thr <- apply(nc, 2L, stats::quantile, 0.95, na.rm = TRUE)
    mem <- which(o$contrib > thr & vapply(obs, length, 1L) >= min_shared)
    if (length(mem) < 2L) next
    v <- setNames(o$v, ids)
    prog <- vapply(cts[mem], function(ct) {
      Z <- fit$scores[[ct]]
      r <- suppressWarnings(abs(stats::cor(Z, v[rownames(Z)]))); colnames(Z)[which.max(r)]
    }, "")
    id <- id + 1L
    out[[id]] <- data.frame(cluster = id, cell_type = cts[mem], program_id = unname(prog), stat = o$value, fwer_p = p,
                            stringsAsFactors = FALSE)
  }
  list(method = "subspace_sum", n_perm = n_perm, clusters = if (length(out)) do.call(rbind, out) else .empty_clusters())
}

## ---- node_agg: node-level evidence aggregation across partner cell types ----
## Port of the cmp_detection finalist `n_agg_noanc` (reports/biomni_transfer/part1_callers_scoring.md
## 1.2; code: benchmarks/cmp_detection/scripts/cand_nodeagg.R). Edge max-T is power-limited when
## cell types share few donors: a program coupled to 3-4 cell types can have every single edge just
## below threshold. Node statistic pools evidence over partner cell types:
##   A_v = sum_{c' != ct(v)} (max_{u in c'} |T_vu| - tau)_+^2,  tau = 2
## on rank-first scores (rank-Gauss, then residualize on strata -- the order matters). Node FWER by
## the permutation max of A_v under the same within-strata row permutations; an edge is kept only
## between two significant nodes with |T| >= link_thr, then greedy single linkage. Under the global
## null P(any call) <= alpha (single-step max over nodes).
.rank_gauss_coord <- function(S) {
  R <- apply(S, 2, function(v) stats::qnorm((rank(v, ties.method = "average") - 0.5) / length(v)))
  if (is.null(dim(R))) R <- matrix(R, nrow(S))
  dimnames(R) <- dimnames(S)
  R
}
.prep_rank_first_coord <- function(fit, strata) {
  lapply(fit$scores, function(S) {
    if (is.null(S) || !ncol(S)) return(S)
    R0 <- .rank_gauss_coord(S)
    f <- factor(strata[rownames(S)])
    X <- if (nlevels(f) > 1L) stats::model.matrix(~f) else matrix(1, nrow(S), 1)
    R <- R0 - X %*% qr.solve(X, R0)
    dimnames(R) <- dimnames(S)
    R
  })
}
method_node_agg <- function(fit, strata, n_perm = 999L, seed = 1L, alpha = 0.05,
                            tau = 2, link_thr = 3, min_shared = 20L, min_density = 0) {
  set.seed(seed)
  Z <- .prep_rank_first_coord(fit, strata)
  Z <- Z[vapply(Z, function(S) !is.null(S) && ncol(S) > 0L, TRUE)]
  nodes <- do.call(rbind, lapply(names(Z), function(ct) {
    if (!ncol(Z[[ct]])) return(NULL)
    data.frame(cell_type = ct, program_id = colnames(Z[[ct]]), stringsAsFactors = FALSE)
  }))
  cts <- names(Z)
  pairs <- list()
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    ids <- intersect(rownames(Z[[cts[i]]]), rownames(Z[[cts[j]]]))
    if (length(ids) >= min_shared) pairs[[length(pairs) + 1L]] <- list(a = cts[i], b = cts[j], ids = ids)
  }
  if (!length(pairs)) return(list(method = "node_agg", n_perm = n_perm, clusters = .empty_clusters()))
  edge_T <- function(p, pm = NULL) {
    Za <- if (is.null(pm)) Z[[p$a]] else pm[[p$a]]; Zb <- if (is.null(pm)) Z[[p$b]] else pm[[p$b]]
    r <- suppressWarnings(stats::cor(Za[p$ids, , drop = FALSE], Zb[p$ids, , drop = FALSE]))
    r[!is.finite(r)] <- 0
    atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(length(p$ids) - 3)
  }
  edges_all <- do.call(rbind, lapply(pairs, function(p) {
    t <- edge_T(p)
    data.frame(ct_a = p$a, ct_b = p$b, prog_a = rownames(t)[row(t)], prog_b = colnames(t)[col(t)],
               stat = as.vector(abs(t)), stringsAsFactors = FALSE)
  }))
  ka <- paste(edges_all$ct_a, edges_all$prog_a, sep = "|")
  kb <- paste(edges_all$ct_b, edges_all$prog_b, sep = "|")
  # A_v per node: max |T| per (node, partner cell type), shift by tau, square, sum over partners
  A_v <- function(stat) {
    key <- paste(c(ka, kb), c(edges_all$ct_b, edges_all$ct_a), sep = "@")
    m <- tapply(c(stat, stat), key, max)
    tapply(pmax(0, m - tau)^2, sub("@.*", "", names(m)), sum)
  }
  A <- A_v(edges_all$stat)
  mx <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    Zp <- lapply(Z, function(S) .perm_rows(S, strata[rownames(S)]))
    stat_p <- unlist(lapply(pairs, function(p) as.vector(abs(edge_T(p, Zp)))), use.names = FALSE)
    mx[b] <- max(A_v(stat_p))
  }
  node_p <- (1 + vapply(A, function(a) sum(mx >= a), numeric(1L))) / (n_perm + 1L)
  sig <- node_p < alpha
  keep <- sig[paste(edges_all$ct_a, edges_all$prog_a, sep = "|")] &
    sig[paste(edges_all$ct_b, edges_all$prog_b, sep = "|")] & edges_all$stat >= link_thr
  edges <- edges_all[keep, , drop = FALSE]
  if (!nrow(edges)) return(list(method = "node_agg", n_perm = n_perm, clusters = .empty_clusters()))
  # cluster fwer_p = max of the two endpoint node p's (part1 1.2 step 6)
  edges$fwer_p <- pmax(node_p[paste(edges$ct_a, edges$prog_a, sep = "|")],
                       node_p[paste(edges$ct_b, edges$prog_b, sep = "|")])
  list(method = "node_agg", n_perm = n_perm, clusters = .cluster_edges(edges, nodes, min_density))
}

METHODS <- c("pa_cc", "node_cor", "node_agg", "subspace_sum", "gene_match")
run_method <- function(name, fit, strata, n_perm = 999L, seed = 1L, ...) {
  switch(name,
         pa_cc = method_pa_cc(fit, strata, n_perm, seed),
         node_cor = method_node_cor(fit, strata, n_perm, seed, ...),
         node_agg = method_node_agg(fit, strata, n_perm, seed, ...),
         subspace_sum = method_subspace_sum(fit, strata, n_perm, seed),
         gene_match = method_gene_match(fit, n_perm, seed),
         stop("unknown method ", name))
}
