## stability.R — subsampling stability of stage-1 programs
##
## Subsamples donors within each cell type (without replacement: duplicated
## rows create spurious low-rank structure for greedy EBMF), refits with the
## stored settings and the same gene sets/covariates, and matches refit
## factors to the originals via match_programs().

#' Assess program stability by donor subsampling
#'
#' @param fit A canonicalized `cell_program_fit` (refit uses its stored
#'   preprocessing settings, gene sets and covariates).
#' @param n_boot Number of subsample refits.
#' @param resample_frac Fraction of donors kept per cell type per refit.
#' @param match_method Matching method for `match_programs`.
#' @param sim_threshold Minimum |cosine| for a factor-level match.
#' @param top_n Top-gene set size for Jaccard recovery.
#' @param seed Random seed.
#' @param progress Print progress messages.
#' @return The fit with `fit$stability` populated: per-program recovery
#'   frequency, median score correlation, median loading cosine, sign
#'   consistency, and top-gene Jaccard.
#' @export
assess_program_stability <- function(fit, n_boot = 100L,
                                     match_method = c("subspace", "correlation"),
                                     sim_threshold = 0.5, top_n = 20L,
                                     resample_frac = 0.8,
                                     seed = NULL, progress = TRUE) {
  match_method <- match.arg(match_method)
  if (!is.null(seed)) set.seed(seed)
  pre <- fit$preprocessing
  prog_ids <- .program_table(fit)$program_id
  rec <- stats::setNames(vector("list", length(prog_ids)), prog_ids)
  for (p in prog_ids) {
    rec[[p]] <- list(hit = 0L, score_cor = c(), loading_cos = c(),
                     sign_consistent = 0L, top_jaccard = c())
  }

  for (b in seq_len(n_boot)) {
    if (progress && b %% 10L == 0L) message(sprintf("subsample %d/%d", b, n_boot))
    x_b <- fit[c("sample_metadata", "feature_metadata", "cell_types", "observation_ids", "shared_genes", "dims")]
    x_b$pseudobulk <- lapply(fit$pseudobulk, function(Y) {
      Y[sample.int(nrow(Y), max(5L, floor(resample_frac * nrow(Y)))), , drop = FALSE]
    })
    class(x_b) <- "cell_program_data"
    fit_b <- tryCatch(
      fit_celltype_programs(
        x_b, loading_prior = pre$loading_prior, max_factors = pre$max_factors,
        center = pre$centered, scale = pre$scaled, features = pre$feature_sets,
        var_type = pre$var_type, backfit = pre$backfit, nullcheck = pre$nullcheck,
        covariates = pre$covariates, seed = fit$provenance$seed + b
      ) |> canonicalize_programs(),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next

    for (ct in names(fit$fits)) {
      K0 <- ncol(fit$scores[[ct]])
      if (K0 == 0L) next
      mp <- match_programs(fit, fit_b, ct, method = match_method,
                           sim_threshold = sim_threshold)
      if (nrow(mp$matches) == 0L) next
      for (i in seq_len(nrow(mp$matches))) {
        pa <- mp$matches$program_a[i]; pb <- mp$matches$program_b[i]
        r <- rec[[pa]]
        r$hit <- r$hit + 1L
        r$loading_cos <- c(r$loading_cos, abs(mp$matches$cosine[i]))
        za <- fit$scores[[ct]][, pa]
        zb <- fit_b$scores[[ct]][, pb]
        common <- intersect(names(za), names(zb))
        if (length(common) >= 5L) {
          r$score_cor <- c(r$score_cor, abs(stats::cor(za[common], zb[common])))
        }
        r$sign_consistent <- r$sign_consistent + as.integer(mp$matches$cosine[i] > 0)
        ta <- top_genes(fit, pa, n = top_n)$gene
        tb <- top_genes(fit_b, pb, n = top_n)$gene
        r$top_jaccard <- c(r$top_jaccard,
                           length(intersect(ta, tb)) / length(union(ta, tb)))
        rec[[pa]] <- r
      }
    }
  }

  stab <- do.call(rbind, lapply(prog_ids, function(p) {
    r <- rec[[p]]
    data.frame(
      program_id = p,
      cell_type = .program_table(fit)$cell_type[match(p, prog_ids)],
      recovery_freq = r$hit / n_boot,
      median_score_cor = if (length(r$score_cor)) stats::median(r$score_cor) else NA_real_,
      median_loading_cos = if (length(r$loading_cos)) stats::median(r$loading_cos) else NA_real_,
      sign_consistency = if (r$hit > 0L) r$sign_consistent / r$hit else NA_real_,
      top_gene_jaccard = if (length(r$top_jaccard)) stats::median(r$top_jaccard) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  fit$stability <- list(per_program = stab, n_boot = n_boot,
                        match_method = match_method, sim_threshold = sim_threshold)
  fit
}
