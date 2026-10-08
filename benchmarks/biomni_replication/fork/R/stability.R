## stability.R — bootstrap stability of stage-1 programs
##
## Bootstrap resamples biological observations within each cell type, refits
## with the stored stage-1 settings (L1/L2 borrowing disabled inside
## bootstrap for speed), and matches refit factors to the original programs
## via subspace-aware matching.

#' Assess program stability by bootstrap
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param n_boot Number of bootstrap refits.
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
  params <- fit$provenance$parameters

  ## per-program collectors
  prog_ids <- fit$programs$program_id
  rec <- stats::setNames(vector("list", length(prog_ids)), prog_ids)
  for (p in prog_ids) {
    rec[[p]] <- list(hit = 0L, score_cor = c(), loading_cos = c(),
                     sign_consistent = 0L, top_jaccard = c())
  }

  for (b in seq_len(n_boot)) {
    if (progress && b %% 10L == 0L) message(sprintf("bootstrap %d/%d", b, n_boot))
    ## subsample observations per cell type (without replacement) and refit.
    ## With-replacement bootstrap duplicates rows, which creates spurious
    ## low-rank structure for greedy EBMF (observed: 10 junk factors vs 2 true,
    ## ~40x slowdown) — subsampling avoids the pathology (Meinshausen-Buhlmann
    ## stability-selection style).
    mats_b <- lapply(fit$input$matrices, function(Y) {
      idx <- sample.int(nrow(Y), size = max(5L, floor(resample_frac * nrow(Y))),
                        replace = FALSE)
      Y[idx, , drop = FALSE]
    })
    x_b <- as_cell_program_data(mats_b)
    fit_b <- tryCatch(
      fit_celltype_programs(
        x_b,
        loading_prior = params$loading_prior, score_prior = params$score_prior,
        max_factors = params$max_factors, center = TRUE, scale = FALSE,
        var_type = params$var_type, backfit = params$backfit,
        nullcheck = params$nullcheck, share_prior = "none",
        verbose = 0L
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
        ## score correlation over shared observations (bootstrap ids differ;
        ## use the original observation prefix)
        za <- fit$scores[[ct]][, pa]
        zb <- fit_b$scores[[ct]][, pb]
        base_ids <- sub("__b\\d+_\\d+$", "", names(zb))
        common <- intersect(names(za), base_ids)
        if (length(common) >= 5L) {
          r$score_cor <- c(r$score_cor,
                           abs(stats::cor(za[common], zb[match(common, base_ids)])))
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
      cell_type = fit$programs$cell_type[match(p, prog_ids)],
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
