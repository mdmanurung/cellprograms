## methods.R -- candidate methods for recurring programs across cell types.
## v2: node_cor permutation fixed (positional, pair-specific; Defect D);
## .maxT_edges perm_fun is now pair-specific; block_E and eb_cp arms added;
## oracle (true-loadings projection) added as the detectability floor.
## Each method returns list(method, clusters = data.frame(cluster, cell_type,
## program_id, stat, fwer_p), activity = donors x clusters matrix or NULL).
## program_id is NA for pa_cc (cell-type-level clusters).

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
.empty_clusters <- function() data.frame(cluster = integer(), cell_type = character(),
                                         program_id = character(),
                                         stat = numeric(), fwer_p = numeric(), stringsAsFactors = FALSE)

## all cross-cell-type node pairs: statistic function applied per cell-type pair, max-T null.
## perm_fun(i) returns the pair-specific permutation payload for pair i.
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
    mx[b] <- max(vapply(seq_along(pairs), function(i) {
      if (is.null(obs[[i]])) return(0)
      max(abs(stat_fun(pairs[[i]], perm_fun(i))))
    }, numeric(1L)))
  }
  edges$fwer_p <- vapply(edges$stat, function(s) (1 + sum(mx >= s)) / (n_perm + 1), numeric(1L))
  edges[edges$fwer_p < alpha, , drop = FALSE]
}

## mean of member node scores per donor (union, available-case), for decoy checks
.cluster_activity <- function(clusters, Z) {
  if (!nrow(clusters)) return(NULL)
  ids <- sort(unique(unlist(lapply(Z, rownames))))
  out <- matrix(NA_real_, length(ids), length(unique(clusters$cluster)),
                dimnames = list(ids, paste0("cluster_", sort(unique(clusters$cluster)))))
  for (cl in sort(unique(clusters$cluster))) {
    mem <- clusters[clusters$cluster == cl, , drop = FALSE]
    acc <- matrix(NA_real_, length(ids), nrow(mem))
    for (k in seq_len(nrow(mem))) {
      z <- Z[[mem$cell_type[k]]]
      if (!mem$program_id[k] %in% colnames(z)) next
      acc[match(rownames(z), ids), k] <- z[, mem$program_id[k]]
    }
    out[, paste0("cluster_", cl)] <- rowMeans(acc, na.rm = TRUE)
  }
  out
}

## ---- node_cor core: studentized score correlation, max-T over node pairs --
## Positional pre-alignment per pair; permutation permutes ROW ORDER of the
## a-side matrix within strata (Defect D fix: never re-index by donor name
## after permuting).
.node_cor_core <- function(Z, strata, n_perm, seed, alpha, min_shared) {
  set.seed(seed)
  cts <- names(Z)
  pairs <- list()
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    a <- cts[i]; b <- cts[j]
    ids <- intersect(rownames(Z[[a]]), rownames(Z[[b]]))
    if (length(ids) < min_shared) next
    pairs[[length(pairs) + 1L]] <- list(
      a = a, b = b,
      Za = Z[[a]][ids, , drop = FALSE],
      Zb = Z[[b]][ids, , drop = FALSE],
      st = unname(strata[ids]))
  }
  stat_fun <- function(p, pm = NULL) {
    Za <- if (is.null(pm)) p$Za else pm
    r <- suppressWarnings(stats::cor(Za, p$Zb))
    r[!is.finite(r)] <- 0
    atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(nrow(p$Za) - 3)
  }
  perm_fun <- function(i) {
    p <- pairs[[i]]
    p$Za[.perm_within(nrow(p$Za), p$st), , drop = FALSE]
  }
  edges <- .maxT_edges(pairs, stat_fun, perm_fun, n_perm, alpha)
  nodes <- do.call(rbind, lapply(cts, function(ct)
    data.frame(cell_type = ct, program_id = colnames(Z[[ct]]), stringsAsFactors = FALSE)))
  cl <- if (is.null(edges) || !nrow(edges)) .empty_clusters() else .cluster_edges(edges, nodes)
  list(method = "node_cor", n_perm = n_perm, clusters = cl,
       activity = .cluster_activity(cl, Z), edges = edges)
}

method_node_cor <- function(fit, strata, n_perm = 999L, seed = 1L, alpha = 0.05, min_shared = 20L) {
  cts <- fit$cell_types[vapply(fit$scores, ncol, 1L) > 0L]
  .node_cor_core(fit$scores[cts], strata, n_perm, seed, alpha, min_shared)
}

## ---- oracle: project adjusted matrices onto TRUE loadings -----------------
## Detectability floor: bypasses stage-1 estimation; nodes are the true
## programs' projection scores. Noise comes only from the data, so power
## increases with amplitude. Runs the same node_cor machinery.
method_oracle <- function(mats, truth, strata, n_perm = 999L, seed = 1L,
                          alpha = 0.05, min_shared = 20L) {
  Z <- list()
  for (ct in names(mats)) {
    cols <- list()
    for (p in truth$programs) if (ct %in% p$members) {
      g <- p$genes[[ct]]; w <- p$w[[ct]]
      cols[[p$id]] <- drop(mats[[ct]][, g, drop = FALSE] %*% w / sum(w^2))
    }
    for (k in seq_along(truth$privates[[ct]])) {
      pv <- truth$privates[[ct]][[k]]
      cols[[paste0(ct, "_P", k)]] <- drop(mats[[ct]][, pv$genes, drop = FALSE] %*% pv$w / sum(pv$w^2))
    }
    if (length(cols)) {
      M <- do.call(cbind, cols)
      rownames(M) <- rownames(mats[[ct]])
      Z[[ct]] <- M
    }
  }
  r <- .node_cor_core(Z, strata, n_perm, seed, alpha, min_shared)
  r$method <- "oracle"
  r
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
    Wb <- if (is.null(pm)) unit(W[[p$b]][p$g, , drop = FALSE]) else unit(W[[p$b]][pm, , drop = FALSE])
    crossprod(Wa, Wb) * sqrt(length(p$g))
  }
  perm_fun <- function(i) match(sample(pairs[[i]]$g), rownames(W[[pairs[[i]]$b]]))
  edges <- .maxT_edges(pairs, stat_fun, perm_fun, n_perm, alpha)
  list(method = "gene_match", n_perm = n_perm, activity = NULL,
       clusters = if (is.null(edges) || !nrow(edges)) .empty_clusters() else .cluster_edges(edges, nodes))
}

## ---- pa_cc: existing pairwise test, connected components over cell types ----
method_pa_cc <- function(fit, strata, n_perm = 1999L, seed = 1L) {
  ss <- sharing_spectrum(fit, space = "scores", n_perm = n_perm, seed = seed, strata = strata)$summary
  called <- ss[!is.na(ss$shared) & ss$shared, , drop = FALSE]
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
  list(method = "pa_cc", n_perm = n_perm, clusters = cl, activity = NULL)
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
  out <- list(); id <- 0L; acts <- list()
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
    acts[[id]] <- v
  }
  cl <- if (length(out)) do.call(rbind, out) else .empty_clusters()
  activity <- NULL
  if (length(acts)) {
    activity <- do.call(cbind, acts)
    colnames(activity) <- paste0("cluster_", seq_along(acts))
  }
  list(method = "subspace_sum", n_perm = n_perm, clusters = cl, activity = activity)
}

## ---- block_E: joint block EBMF at expression level -------------------------
method_block_e <- function(fit, max_factors = 20L, seed = 1L, mass_thresh = 0.1,
                           strata = NULL, n_perm = 0L, alpha = 0.05) {
  bfit <- fit_joint_block(fit, level = "expression", order = "shared_first",
                          max_factors = max_factors, seed = seed, verbose = 0L)
  bc <- block_clusters(bfit, mass_thresh = mass_thresh)
  out <- list(method = "block_E", n_perm = n_perm, clusters = bc$clusters,
              activity = bc$activity, fit = bfit)
  ## Calibration (post-prereg development iteration, decision-table follow-up):
  ## the block-mass membership rule has no significance layer -- on null data the
  ## shared-first pass still returns ~Kmax factors whose loadings spread across
  ## blocks (main run: P(any spurious) = 1.00, mean 15.9 calls/null). With
  ## n_perm > 0 and strata supplied, keep only factors whose donor-level signal
  ## is significant under within-strata donor permutation applied independently
  ## per cell-type block (destroys cross-block coordination, preserves
  ## within-block structure); max-T over factors -> FWER.
  if (n_perm > 0L && !is.null(strata) && bfit$shared_fit$n_factors > 0L &&
      nrow(bc$clusters)) {
    ## Split-half calibration (eb_cp-style held-out selection): the full-data
    ## fit's factors are optimized on the same donors the permutation tests, so
    ## noise factors inherit winner's-curse alignment and survive a fixed-
    ## loading null (probe round 3: null FPR 7/18, true call lost). Instead:
    ## refit the joint pass on half A (stratified), test cross-block alignment
    ## of the half-A factors on held-out half-B donors (no winner's curse by
    ## construction), and report the half-A fit's calibrated calls.
    donors <- sort(unique(unlist(lapply(fit$fits[fit$cell_types],
                                        function(f) rownames(f$Y_used)))))
    set.seed(seed)
    inA <- unlist(lapply(split(donors, strata[donors]), function(g)
      sample(g, ceiling(length(g) / 2))))
    inA <- sort(inA)
    fitA <- fit
    fitA$fits <- lapply(fit$fits, function(f) {
      f$Y_used <- f$Y_used[rownames(f$Y_used) %in% inA, , drop = FALSE]; f
    })
    bfitA <- fit_joint_block(fitA, level = "expression", order = "shared_first",
                             max_factors = max_factors, seed = seed, verbose = 0L)
    bcA <- block_clusters(bfitA, mass_thresh = mass_thresh)
    if (bfitA$shared_fit$n_factors > 0L && nrow(bcA$clusters)) {
      pm <- .block_mass_maxT(bfitA, fit, strata, n_perm, seed,
                             mass_thresh = mass_thresh,
                             heldout = setdiff(donors, inA))
      out$factor_p <- pm
      keep_id <- pm$factor_id[pm$p < alpha]
      keep <- bcA$clusters$program_id %in% keep_id
      out$clusters <- bcA$clusters[keep, , drop = FALSE]
      ## activity: project FULL data onto kept factors, donor-mean over members
      kept_cl <- unique(bcA$clusters$cluster[keep])
      if (length(kept_cl)) {
        fids <- sprintf("SF%02d", seq_len(bfitA$shared_fit$n_factors))
        acts <- list()
        for (cid in kept_cl) {
          pid <- bcA$clusters$program_id[bcA$clusters$cluster == cid][1]
          k <- match(pid, fids)
          mem <- unique(bcA$clusters$cell_type[bcA$clusters$cluster == cid])
          vs <- lapply(mem, function(ct) {
            sl <- bfitA$shared_fit$F_pm[bfitA$block_id == ct, k]
            setNames(as.vector(fit$fits[[ct]]$Y_used %*% sl),
                     rownames(fit$fits[[ct]]$Y_used))
          })
          all_d <- sort(unique(unlist(lapply(vs, names))))
          acts[[paste0("cluster_", cid)]] <- vapply(all_d, function(d)
            mean(vapply(vs, function(v) if (d %in% names(v)) v[[d]] else NA_real_,
                        numeric(1)), na.rm = TRUE), numeric(1))
        }
        activity <- do.call(cbind, acts)
        colnames(activity) <- names(acts)
        out$activity <- activity
      }
    } else {
      out$clusters <- bcA$clusters[0, , drop = FALSE]
      out$activity <- NULL
    }
  }
  out
}

## Permutation max-T p-values for block_E shared factors. Statistic for factor
## k: the mean pairwise |cor| of the factor's donor activities across the block
## pairs that pass the block-mass rule for k (the pairs a shared call would
## rest on); factors with < 2 mass-passing blocks cannot produce a call and get
## p = 1. Under the null of no cross-block coordination, each block's donor rows
## are permuted within strata independently, which destroys cross-block activity
## alignment while leaving within-block (private) structure intact. A
## total-projection statistic ||X f_k||^2 is NOT usable here: it is exactly
## invariant to within-strata row permutation for single-block factors, so
## every private factor would be "significant" by construction (observed in the
## first probe round: most factors pinned at the p = 1/(n_perm+1) floor). f_k is
## fixed from the observed fit (winner's-curse caveat: f_k was selected to
## explain variance on this data; the cross-block alignment statistic is not
## directly optimized by the fit, and the probes measure the realized FPR).
.block_mass_maxT <- function(bfit, fit, strata, n_perm = 199L, seed = 1L,
                            mass_thresh = 0.1, heldout = NULL) {
  F_pm <- bfit$shared_fit$F_pm                    # stacked genes x K
  K <- ncol(F_pm)
  fids <- sprintf("SF%02d", seq_len(K))           # matches fit_joint_block summary ids
  cts <- colnames(bfit$block_mass)
  blocks <- lapply(fit$fits[cts], function(f) {
    Y <- f$Y_used
    if (!is.null(heldout)) Y <- Y[rownames(Y) %in% heldout, , drop = FALSE]
    Y
  })
  keep_ct <- vapply(blocks, function(Y) nrow(Y) >= 20L, TRUE)
  blocks <- blocks[keep_ct]
  cts <- cts[keep_ct]
  mass <- bfit$block_mass[, cts, drop = FALSE]
  ## loading slice of every factor on each block's genes (assert alignment)
  slices <- lapply(setNames(cts, cts), function(ct) {
    S <- F_pm[bfit$block_id == ct, , drop = FALSE]
    rn <- sub(paste0("^", ct, "__"), "", rownames(S))
    stopifnot(identical(rn, colnames(blocks[[ct]])))
    S
  })
  sts <- lapply(blocks, function(Y) factor(strata[rownames(Y)]))
  mass <- bfit$block_mass                         # K x cts
  pairs <- combn(cts, 2L)
  ## donor activities per block and factor; rows matched by donor id
  vof <- function(bl) lapply(setNames(cts, cts), function(ct) bl[[ct]] %*% slices[[ct]])
  ## T_k: mean |cor| across mass-passing block pairs (0 if < 2 such blocks)
  Tof <- function(vs) {
    vapply(seq_len(K), function(k) {
      mem <- cts[mass[k, ] > mass_thresh]
      if (length(mem) < 2L) return(0)
      pr <- pairs[, pairs[1, ] %in% mem & pairs[2, ] %in% mem, drop = FALSE]
      if (ncol(pr) == 0L) return(0)
      mean(vapply(seq_len(ncol(pr)), function(j) {
        a <- vs[[pr[1, j]]]; b <- vs[[pr[2, j]]]
        common <- intersect(rownames(a), rownames(b))
        if (length(common) < 10L) return(0)
        x <- a[common, k]; y <- b[common, k]
        if (stats::sd(x) < 1e-12 || stats::sd(y) < 1e-12) 0 else suppressWarnings(cor(x, y))
      }, numeric(1)))
    }, numeric(1))
  }
  T_obs <- Tof(vof(blocks))
  set.seed(seed)
  mx <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    bl_p <- lapply(setNames(cts, cts), function(ct) {
      Y <- blocks[[ct]]
      Yp <- Y[.perm_within(nrow(Y), sts[[ct]]), , drop = FALSE]
      rownames(Yp) <- rownames(Y)   # keep donor labels fixed: name-based donor
                                    # matching must pair PERMUTED values, else the
                                    # rowname-carried permutation is undone
      Yp
    })
    tp <- Tof(vof(bl_p))
    mx[b] <- max(tp[tp > 0], 0)                   # max-T over callable factors
  }
  p <- vapply(T_obs, function(t) if (t <= 0) 1 else (1 + sum(mx >= t)) / (n_perm + 1), numeric(1))
  data.frame(factor_id = fids, T_obs = T_obs, p = p, stringsAsFactors = FALSE)
}

## ---- eb_cp: EB-CP tensor factorization ------------------------------------
method_eb_cp <- function(fit, strata = NULL, Kmax = 15L, seed = 1L, member_thresh = 0.2) {
  mats <- lapply(fit$fits[fit$cell_types], function(f) f$Y_used)
  cp <- eb_cp_fit(mats, Kmax = Kmax, seed = seed, member_thresh = member_thresh,
                  strata = strata)
  cc <- eb_cp_clusters(cp, member_thresh)
  list(method = "eb_cp", n_perm = 0L, clusters = cc$clusters, activity = cc$activity,
       fit = cp)
}

METHODS <- c("pa_cc", "node_cor", "subspace_sum", "gene_match", "block_E", "eb_cp")
run_method <- function(name, fit, strata, n_perm = 999L, seed = 1L) {
  switch(name,
         pa_cc = method_pa_cc(fit, strata, n_perm = 1999L, seed = seed),
         node_cor = method_node_cor(fit, strata, n_perm, seed),
         subspace_sum = method_subspace_sum(fit, strata, n_perm, seed),
         gene_match = method_gene_match(fit, n_perm, seed),
         block_E = method_block_e(fit, seed = seed),
         eb_cp = method_eb_cp(fit, strata = strata, seed = seed),
         stop("unknown method ", name))
}
