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

## permute the rows of M within strata (rows are donors of one cell type)
.perm_rows <- function(M, strata_ct) M[.perm_within(nrow(M), strata_ct), , drop = FALSE]

## Greedy single-linkage on called edges, strongest first; a merge that would put two
## nodes of the same cell type in one cluster is skipped.
.cluster_edges <- function(edges, nodes) {
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
    id <- id + 1L
    e <- edges[used[[nm]], , drop = FALSE]
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
method_node_cor <- function(fit, strata, n_perm = 999L, seed = 1L, alpha = 0.05, min_shared = 20L) {
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
       clusters = if (is.null(edges) || !nrow(edges)) .empty_clusters() else .cluster_edges(edges, nodes))
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
  spec <- function(Us) {
    E <- Map(embed, Us, cts)
    M <- Reduce(`+`, lapply(E, function(e) tcrossprod(e)))
    ev <- eigen(Dh * M * rep(Dh, each = nu), symmetric = TRUE)
    V <- ev$vectors[, seq_len(J), drop = FALSE] * Dh        # back to donor space
    contrib <- vapply(seq_along(cts), function(k) {
      o <- obs[[k]]; v <- V[o, , drop = FALSE]
      colSums((crossprod(E[[k]][o, , drop = FALSE], v))^2) / pmax(colSums(v^2), 1e-12)
    }, numeric(J))
    list(values = ev$values[seq_len(J)], V = V, contrib = contrib)   # contrib: J x n_ct
  }
  o <- spec(Ub)
  nl <- matrix(NA_real_, n_perm, J); nc <- array(NA_real_, c(n_perm, J, length(cts)))
  st <- strata
  for (b in seq_len(n_perm)) {
    Up <- Map(function(U, ct) .perm_rows(U, st[rownames(fit$scores[[ct]])]), Ub, cts)
    s <- spec(Up); nl[b, ] <- s$values; nc[b, , ] <- s$contrib
  }
  p <- vapply(seq_len(J), function(j) (1 + sum(nl[, j] >= o$values[j])) / (n_perm + 1), numeric(1L))
  p <- cummax(p)                                   # step-down: monotone in j
  out <- list(); id <- 0L
  for (j in seq_len(J)) {
    if (p[j] >= alpha) break
    thr <- apply(nc[, j, ], 2L, stats::quantile, 0.95, na.rm = TRUE)
    mem <- which(o$contrib[j, ] > thr & vapply(obs, length, 1L) >= min_shared)
    if (length(mem) < 2L) next
    v <- setNames(o$V[, j], ids)
    prog <- vapply(cts[mem], function(ct) {
      Z <- fit$scores[[ct]]
      r <- suppressWarnings(abs(stats::cor(Z, v[rownames(Z)]))); colnames(Z)[which.max(r)]
    }, "")
    id <- id + 1L
    out[[id]] <- data.frame(cluster = id, cell_type = cts[mem], program_id = unname(prog), stat = o$values[j], fwer_p = p[j],
                            stringsAsFactors = FALSE)
  }
  list(method = "subspace_sum", n_perm = n_perm, clusters = if (length(out)) do.call(rbind, out) else .empty_clusters())
}

METHODS <- c("pa_cc", "node_cor", "subspace_sum", "gene_match")
run_method <- function(name, fit, strata, n_perm = 999L, seed = 1L) {
  switch(name,
         pa_cc = method_pa_cc(fit, strata, n_perm, seed),
         node_cor = method_node_cor(fit, strata, n_perm, seed),
         subspace_sum = method_subspace_sum(fit, strata, n_perm, seed),
         gene_match = method_gene_match(fit, n_perm, seed),
         stop("unknown method ", name))
}
