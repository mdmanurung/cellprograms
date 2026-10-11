## cand_flash.R -- F5: EBMF (flashier) on the donor x node score matrix (rank-gauss, site-residualized,
## NA for donors a cell type lacks). Donor side normal prior, node side point-Laplace (sparse membership).
## A factor is a CMP if >= 2 cell types have a node with lfsr < lfsr_thr; one node per cell type (largest |F|).
## Optional permutation calibration (perm > 0): the factor statistic (sum over cts of top (1-lfsr) F^2 / ||F||^2
## restricted to significant nodes) must exceed the 1-alpha quantile of the max statistic on within-site permuted data.
.node_matrix <- function(Z) {
  ids <- sort(unique(unlist(lapply(Z, rownames))))
  M <- do.call(cbind, lapply(names(Z), function(c) {
    S <- matrix(NA_real_, length(ids), ncol(Z[[c]]), dimnames = list(ids, paste(c, colnames(Z[[c]]), sep = "|")))
    S[rownames(Z[[c]]), ] <- Z[[c]]; S
  }))
  M
}
.flash_cmp <- function(M, Kmax, lfsr_thr) {
  fl <- suppressWarnings(suppressMessages(flashier::flash(M, ebnm_fn = c(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
                                                           greedy_Kmax = Kmax, var_type = 2, backfit = TRUE, verbose = 0)))
  if (is.null(fl$n_factors) || fl$n_factors == 0) return(list())
  F <- fl$F_pm; L <- fl$F_lfsr; if (is.null(L)) L <- matrix(0, nrow(F), ncol(F))
  ct <- sub("\\|.*", "", colnames(M))
  out <- list()
  for (k in seq_len(ncol(F))) {
    sig <- which(L[, k] < lfsr_thr)
    if (!length(sig)) next
    sel <- sig[order(-abs(F[sig, k]))]; sel <- sel[!duplicated(ct[sel])]
    if (length(sel) < 2L) next
    stat <- sum((1 - L[sel, k]) * F[sel, k]^2) / sum(F[, k]^2)
    out[[length(out) + 1L]] <- list(nodes = colnames(M)[sel], stat = stat)
  }
  out
}
cand_flash_nodes <- function(fit, strata, seed = 1L, Kmax = 8L, lfsr_thr = 0.01, perm = 0L, alpha = 0.05, resid_site = TRUE) {
  set.seed(seed)
  Z <- lapply(lapply(.prep_scores(fit, strata, resid_site), .rank_gauss), .std_cols)
  if (length(Z) < 2L) return(.empty_clusters())
  res <- .flash_cmp(.node_matrix(Z), Kmax, lfsr_thr)
  if (!length(res)) return(.empty_clusters())
  thr <- -Inf
  if (perm > 0) {
    mx <- vapply(seq_len(perm), function(b) {
      r <- .flash_cmp(.node_matrix(.perm_Z(Z, strata)), Kmax, lfsr_thr)
      if (length(r)) max(vapply(r, `[[`, 1, "stat")) else 0
    }, 1)
    thr <- stats::quantile(mx, 1 - alpha)
  }
  out <- list(); id <- 0L
  for (r in res) if (r$stat > thr) {
    id <- id + 1L
    out[[id]] <- data.frame(cluster = id, cell_type = sub("\\|.*", "", r$nodes), program_id = sub("^[^|]*\\|", "", r$nodes),
                            stat = r$stat, fwer_p = NA_real_, stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else .empty_clusters()
}
