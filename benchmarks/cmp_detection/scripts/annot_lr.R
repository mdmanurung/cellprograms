## annot_lr.R -- Tier-2 LR annotation of a called CMP ("LR-compatible", not causal).
## For member nodes a (sender) != b (receiver) of different cell types and LR pair (L, R):
##   s(L,R,a,b) = w_a[L] * w_b[R], with w = |loading| * (1 - lfsr) scaled to max 1 within the node
## (both directions considered). Null: degree-preserving resource permutation -- ligand labels are
## permuted among ligands and receptor labels among receptors (same multiset of pair degrees), the
## max score over the CMP is recomputed; per-pair p = P(null max >= s). Pairs with p < alpha are reported.
.node_w <- function(fit, ct, prog) {
  W <- fit$loadings[[ct]]; if (is.null(W) || !prog %in% colnames(W)) return(NULL)
  w <- abs(W[, prog]); lf <- fit$fits[[ct]]$flash$F_lfsr
  if (!is.null(lf)) w <- w * (1 - lf[, match(prog, colnames(W))])
  w / max(w, 1e-12)
}
annotate_lr <- function(fit, cmp, lr, n_perm = 499L, alpha = 0.05, seed = 1L, top = 10L) {
  set.seed(seed)
  if (nrow(cmp) < 2L) return(NULL)
  ws <- lapply(seq_len(nrow(cmp)), function(i) .node_w(fit, cmp$cell_type[i], cmp$program_id[i]))
  names(ws) <- paste(cmp$cell_type, cmp$program_id, sep = "|")
  ws <- ws[!vapply(ws, is.null, TRUE)]
  if (length(ws) < 2L) return(NULL)
  lr <- unique(lr[, c("ligand", "receptor")])
  score_all <- function(lig, rec) {
    best <- list()
    for (a in names(ws)) for (b in names(ws)) {
      if (sub("\\|.*", "", a) == sub("\\|.*", "", b)) next
      sl <- ws[[a]][lig]; sr <- ws[[b]][rec]; s <- sl * sr; s[is.na(s)] <- 0
      best[[length(best) + 1L]] <- data.frame(k = seq_along(lig), sender = a, receiver = b, s = s)
    }
    do.call(rbind, best)
  }
  obs <- score_all(lr$ligand, lr$receptor)
  obs <- obs[obs$s > 0, , drop = FALSE]
  if (!nrow(obs)) return(data.frame())
  ul <- unique(lr$ligand); ur <- unique(lr$receptor)
  mx <- vapply(seq_len(n_perm), function(b) {
    pl <- setNames(sample(ul), ul); pr <- setNames(sample(ur), ur)
    s <- score_all(unname(pl[lr$ligand]), unname(pr[lr$receptor]))$s
    if (length(s)) max(s) else 0
  }, 1)
  obs$p <- vapply(obs$s, function(s) (1 + sum(mx >= s)) / (n_perm + 1), 1)
  obs$ligand <- lr$ligand[obs$k]; obs$receptor <- lr$receptor[obs$k]
  obs <- obs[order(obs$p, -obs$s), c("sender", "receiver", "ligand", "receptor", "s", "p")]
  utils::head(obs[obs$p < alpha, , drop = FALSE], top)
}
