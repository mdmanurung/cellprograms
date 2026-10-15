## cand_ext.R -- F4 variant: consensus extension of edge-called CMPs.
## Stage 1: edge max-T CMPs (rank-first, ambient filter). Stage 2: for each called CMP k form a consensus donor score
## c_k = mean of member node scores (available-case). Every unassigned node v in a cell type not yet in CMP k is tested
## against c_k with T = atanh(r) sqrt(n - 3). Null: permute v's rows within strata with c_k held fixed; max over all
## (k, v) tests per permutation (max-T, FWER across extensions). v's permutation null is valid conditional on c_k because
## c_k is built from other cell types only. Best significant node per (k, cell type) is added; one pass.
.consensus <- function(Z, cl) {
  ids <- sort(unique(unlist(lapply(seq_len(nrow(cl)), function(i) rownames(Z[[cl$cell_type[i]]])))))
  M <- matrix(NA_real_, length(ids), nrow(cl), dimnames = list(ids, NULL))
  for (i in seq_len(nrow(cl))) { S <- Z[[cl$cell_type[i]]]; M[rownames(S), i] <- S[, cl$program_id[i]] }
  c0 <- rowMeans(M, na.rm = TRUE); c0[rowSums(!is.na(M)) == 0] <- NA; c0
}

cand_extend <- function(fit, strata, seed = 1L, n_perm = 199L, alpha = 0.05, alpha_ext = 0.05, amb_thr = 0.3,
                        min_shared = 20L, rank_first = TRUE, stage1 = NULL) {
  cl <- if (is.null(stage1)) cand_edge(fit, strata, seed, stat = "rankgauss", n_perm = n_perm, alpha = alpha,
                                       amb_thr = amb_thr, rank_first = rank_first) else stage1
  if (!nrow(cl)) return(cl)
  set.seed(seed + 1L)
  Z <- if (rank_first) .prep_rank_first(fit, strata, TRUE) else .transform_scores(.prep_scores(fit, strata, TRUE), "rankgauss")
  if (is.finite(amb_thr) && !is.null(fit$profiles)) {
    a <- ambient_scores(fit, fit$profiles); bad <- a[a$amb > amb_thr, ]
    for (i in seq_len(nrow(bad))) if (!is.null(Z[[bad$cell_type[i]]])) Z[[bad$cell_type[i]]] <- Z[[bad$cell_type[i]]][, colnames(Z[[bad$cell_type[i]]]) != bad$program_id[i], drop = FALSE]
    Z <- Z[vapply(Z, ncol, 1L) > 0]
  }
  used <- paste(cl$cell_type, cl$program_id)
  tests <- list()
  for (k in unique(cl$cluster)) {
    ck <- cl[cl$cluster == k, ]; c0 <- .consensus(Z, ck)
    for (ct in setdiff(names(Z), ck$cell_type)) {
      S <- Z[[ct]]; free <- colnames(S)[!paste(ct, colnames(S)) %in% used]
      if (!length(free)) next
      d <- intersect(rownames(S), names(c0)[!is.na(c0)])
      if (length(d) < min_shared) next
      tests[[length(tests) + 1L]] <- list(k = k, ct = ct, progs = free, d = d, c = c0[d])
    }
  }
  if (!length(tests)) return(cl)
  tstat <- function(S, t) { r <- suppressWarnings(stats::cor(S[t$d, t$progs, drop = FALSE], t$c)); r[!is.finite(r)] <- 0
    atanh(pmin(pmax(r, -0.999999), 0.999999))[, 1] * sqrt(length(t$d) - 3) }
  obs <- lapply(tests, function(t) abs(tstat(Z[[t$ct]], t)))
  mx <- vapply(seq_len(n_perm), function(b) {
    Zp <- .perm_Z(Z, strata)
    max(vapply(tests, function(t) max(abs(tstat(Zp[[t$ct]], t))), 1))
  }, 1)
  add <- list()
  for (i in seq_along(tests)) {
    p <- vapply(obs[[i]], function(s) (1 + sum(mx >= s)) / (n_perm + 1), 1)
    j <- which.max(obs[[i]])
    if (p[j] < alpha_ext) add[[length(add) + 1L]] <- data.frame(cluster = tests[[i]]$k, cell_type = tests[[i]]$ct,
      program_id = tests[[i]]$progs[j], stat = obs[[i]][j], fwer_p = p[j], stringsAsFactors = FALSE)
  }
  if (!length(add)) return(cl)
  add <- do.call(rbind, add)
  add <- add[order(-add$stat), ]; add <- add[!duplicated(paste(add$cell_type, add$program_id)), ]   # a node joins one CMP
  rbind(cl, add[, names(cl)])
}
