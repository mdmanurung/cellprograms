## cand_nodeagg_wt.R -- WS2.1 score-uncertainty arm: propagate flashier score posterior SDs
## (L_psd) into the caller by empirical-Bayes shrinkage of donor scores toward stratum means,
## then run the EXACT winner pipeline (cand_node_agg with anchors = FALSE, i.e. n_agg_noanc)
## on the shrunk scores. The comparison against n_agg_noanc therefore isolates uncertainty
## propagation: identical prep, edge statistic, node statistic, permutation null, grouping.
##
## Shrinkage model (per cell type, per program, per stratum): observed canonicalized score
##   S_ik = mu_gk + signal_ik + noise_ik,  noise_ik ~ (0, sd_ik^2),  sd_ik = L_psd[,k] / s_k
## (s_k = canonicalization scale multiplier: canonicalized score = L_pm/s_k*sign, so the SD
## scales by 1/s_k; sign flips leave SDs unchanged). Method-of-moments signal variance
## V_gk = max(var_g(S[,k]) - mean_g(sd_ik^2), floor); EB posterior mean
##   S*_ik = mu_gk + (S_ik - mu_gk) * V_gk / (V_gk + sd_ik^2).
## High-variance donor scores are shrunk to the stratum mean and contribute ~nothing to the
## edge correlations -- the inverse-variance down-weighting of WS2.1. Deterministic prep, so
## the winner's within-strata permutation null remains valid exactly as for the ranks.
##
## Scope: MARGINAL mode only. In within-context mode (fit$scores_raw set by the driver) the
## score/SD alignment is ambiguous (within_context transforms the fit but scores_raw holds
## marginal scores), so the arm falls back to the unmodified winner there.
## Needs cand_edge.R (.prep_rank_first via cand_node_agg), cand_nodeagg.R (cand_node_agg).

# Canonicalized score SDs for one cell type, aligned to the score matrix S (obs x K).
# NULL when unavailable/mismatched (caller then leaves that cell type's scores unchanged).
.score_sd_canon <- function(fit, ct, S) {
  fl <- fit$fits[[ct]]$flash
  if (is.null(fl) || is.null(fl$L_psd)) return(NULL)
  SD <- fl$L_psd
  if (nrow(SD) != nrow(S) || ncol(SD) != ncol(S)) return(NULL)
  if (!is.null(rownames(SD)) && !identical(rownames(SD), rownames(S))) {
    idx <- match(rownames(S), rownames(SD))
    if (anyNA(idx)) return(NULL)
    SD <- SD[idx, , drop = FALSE]
  }
  s_k <- fit$canonicalization[[ct]]$scale_multiplier
  if (is.null(s_k) || length(s_k) != ncol(S) || any(!is.finite(s_k)) || any(s_k <= 0)) return(NULL)
  SD <- sweep(SD, 2L, s_k, "/")
  colnames(SD) <- colnames(S)
  SD
}

# EB shrinkage of score matrix S (obs x K) toward stratum means, per program.
.eb_shrink_scores <- function(S, SD, strata_ct, min_stratum = 5L, var_floor = 1e-8) {
  f <- factor(strata_ct)
  out <- S
  shrink_block <- function(idx) {
    if (length(idx) < min_stratum) idx <- seq_len(nrow(S))   # small stratum -> global moments
    m <- colMeans(S[idx, , drop = FALSE])
    noise2 <- colMeans(SD[idx, , drop = FALSE]^2)
    sig2 <- pmax(apply(S[idx, , drop = FALSE], 2L, stats::var) - noise2, var_floor)
    rel <- outer(rep(1, nrow(S)), sig2) / (outer(rep(1, nrow(S)), sig2) + SD^2)
    outer(rep(1, nrow(S)), m) + (S - outer(rep(1, nrow(S)), m)) * rel
  }
  for (g in levels(f)) {
    idx <- which(f == g)
    if (!length(idx)) next
    out[idx, ] <- shrink_block(idx)[idx, , drop = FALSE]
  }
  dimnames(out) <- dimnames(S)
  out
}

cand_nodeagg_eb <- function(fit, strata, seed = 1L, tau = 2, link_thr = 3, n_perm = 199L,
                            alpha = 0.05, min_shared = 20L, amb_thr = 0.3) {
  if (is.null(fit$scores_raw)) {                 # marginal mode: shrink canonicalized scores
    shrunk <- fit$scores
    for (ct in names(shrunk)) {
      S <- shrunk[[ct]]
      if (is.null(S) || ncol(S) == 0L) next
      SD <- .score_sd_canon(fit, ct, S)
      if (is.null(SD)) next                      # no usable SDs -> winner behaviour for this ct
      shrunk[[ct]] <- .eb_shrink_scores(S, SD, strata[rownames(S)])
    }
    fit$scores <- shrunk
  }                                              # else within mode: leave scores_raw untouched
  cand_node_agg(fit, strata, seed, tau = tau, link_thr = link_thr, n_perm = n_perm,
                alpha = alpha, min_shared = min_shared, amb_thr = amb_thr, anchors = FALSE)
}
