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
#' @return The fit with `fit$stability` populated:
#'   \describe{
#'     \item{per_program}{One row per program: `n_hit` (refits whose match
#'       recovered the program), `n_eligible` (successful refits in which the
#'       program's cell type produced at least one factor, i.e. matching was
#'       possible), `recovery_freq` (**legacy**: unconditional hit fraction
#'       over all REQUESTED iterations, `n_hit / n_boot` — not an estimate of
#'       biological recovery conditional on successful refitting),
#'       `conditional_recovery_freq` (`n_hit / valid_refits`, NA when no
#'       refit succeeded), median score correlation, median loading cosine,
#'       sign consistency, and top-gene Jaccard.}
#'     \item{valid_refits, failed_refits}{Refit denominator accounting.
#'       Failed refits are never silently dropped: any failure warns, and an
#'       all-fail run yields NA rates (undefined), never an artificial 0.}
#'     \item{replicate_log}{One row per requested refit: `replicate`, `seed`,
#'       `status` ("ok"/"error"), `error_class`/`error_message`,
#'       `elapsed_seconds`, `n_obs` (sampled observations per cell type),
#'       `n_factors` (refit factor count per cell type), and `eligible`
#'       (whether matching could be attempted). A successful refit with no
#'       matched factor is distinguishable from a failed refit.}
#'   }
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' # few subsamples for a quick example; use more (e.g. 100) in practice
#' fit <- assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE)
#' fit$stability$per_program
#' @export
assess_program_stability <- function(fit, n_boot = 100L,
                                     match_method = c("subspace", "correlation"),
                                     sim_threshold = 0.5, top_n = 20L,
                                     resample_frac = 0.8,
                                     seed = NULL, progress = TRUE) {
  match_method <- match.arg(match_method)
  if (!is.null(seed)) set.seed(seed)
  pre <- fit$preprocessing
  prog_tab <- .program_table(fit)
  prog_ids <- prog_tab$program_id
  rec <- stats::setNames(vector("list", length(prog_ids)), prog_ids)
  for (p in prog_ids) {
    rec[[p]] <- list(hit = 0L, eligible = 0L, score_cor = c(), loading_cos = c(),
                     sign_consistent = 0L, top_jaccard = c())
  }

  failed_refits <- 0L
  repl_log <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    if (progress && b %% 10L == 0L) message(sprintf("subsample %d/%d", b, n_boot))
    t0 <- Sys.time()
    # Subsample rows per cell type, then rebuild a CONSISTENT data object:
    # observation_ids, sample_metadata and dims are subset to the sampled
    # rows (previously the full-data values were retained, leaving dims and
    # metadata stale relative to the subsampled matrices).
    sampled_ids <- lapply(fit$pseudobulk, function(Y) {
      Y[sample.int(nrow(Y), max(5L, floor(resample_frac * nrow(Y)))), , drop = FALSE]
    })
    keep <- unique(unlist(lapply(sampled_ids, rownames)))
    x_b <- list(
      pseudobulk = sampled_ids,
      sample_metadata = fit$sample_metadata[fit$sample_metadata$observation_id %in% keep, , drop = FALSE],
      feature_metadata = fit$feature_metadata,
      cell_types = fit$cell_types,
      observation_ids = keep,
      shared_genes = fit$shared_genes,
      dims = lapply(sampled_ids, dim)
    )
    class(x_b) <- "cell_program_data"
    n_obs <- vapply(sampled_ids, nrow, integer(1))
    # CAND-12: replay the frozen preprocessing spec through the MR-0
    # fail-closed checks when the fit carries one, instead of re-deriving
    # preprocessing from recorded settings (which bypassed the checks).
    # Legacy fits without a stored spec keep the explicit-args path, now with
    # feature_scale passed through (it was silently "unknown" in refits).
    spec <- fit$preprocessing$spec
    fit_args <- if (is.null(spec)) {
      list(x_b, loading_prior = pre$loading_prior, max_factors = pre$max_factors,
           center = pre$centered, scale = pre$scaled, features = pre$feature_sets,
           var_type = pre$var_type, backfit = pre$backfit, nullcheck = pre$nullcheck,
           covariates = pre$covariates, covariate_mode = pre$covariate_mode,
           flash_control = pre$flash_control, feature_scale = pre$feature_scale,
           seed = fit$provenance$seed + b)
    } else {
      list(x_b, preprocessing_spec = spec,
           loading_prior = pre$loading_prior, max_factors = pre$max_factors,
           var_type = pre$var_type, backfit = pre$backfit, nullcheck = pre$nullcheck,
           flash_control = pre$flash_control,
           backend = if (is.null(pre$backend)) "flashier" else pre$backend,
           var_floor = if (is.null(pre$var_floor)) "none" else pre$var_floor,
           seed = fit$provenance$seed + b)
    }
    # CP-1: failures are captured with their condition class/message into the
    # replicate log — never collapsed invisibly into zero scientific recovery.
    refit_err <- NULL
    fit_b <- tryCatch(
      do.call(fit_celltype_programs, fit_args) |> canonicalize_programs(),
      error = function(e) { refit_err <<- e; NULL }
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(fit_b)) {
      failed_refits <- failed_refits + 1L
      repl_log[[b]] <- data.frame(
        replicate = b, seed = fit$provenance$seed + b, status = "error",
        error_class = class(refit_err)[1],
        error_message = substr(conditionMessage(refit_err), 1, 200),
        elapsed_seconds = elapsed,
        n_obs = I(list(n_obs)),
        n_factors = I(list(stats::setNames(integer(0), character(0)))),
        eligible = FALSE, stringsAsFactors = FALSE
      )
      next
    }
    k_by_ct <- vapply(fit_b$scores, function(s) if (is.null(s)) 0L else ncol(s), integer(1))
    repl_log[[b]] <- data.frame(
      replicate = b, seed = fit$provenance$seed + b, status = "ok",
      error_class = NA_character_, error_message = NA_character_,
      elapsed_seconds = elapsed,
      n_obs = I(list(n_obs)), n_factors = I(list(k_by_ct)),
      eligible = any(k_by_ct > 0L), stringsAsFactors = FALSE
    )

    for (ct in names(fit$fits)) {
      K0 <- ncol(fit$scores[[ct]])
      if (K0 == 0L) next
      ct_eligible <- isTRUE(k_by_ct[[ct]] > 0L)
      mp <- NULL
      if (ct_eligible) {
        mp <- match_programs(fit, fit_b, ct, method = match_method,
                             sim_threshold = sim_threshold)
      }
      for (p in prog_ids[prog_tab$cell_type == ct]) {
        if (ct_eligible) rec[[p]]$eligible <- rec[[p]]$eligible + 1L
      }
      if (is.null(mp) || nrow(mp$matches) == 0L) next
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

  # CP-1 denominator contract: recovery_freq keeps its legacy meaning
  # (unconditional hit fraction over REQUESTED iterations);
  # conditional_recovery_freq is defined over SUCCESSFUL refits. An all-fail
  # run yields NA for both rates (undefined), never an artificial 0 that
  # would read as "zero stability".
  valid_refits <- as.integer(n_boot) - failed_refits
  if (failed_refits == n_boot) {
    warning(sprintf(paste0("assess_program_stability(): all %d subsample refits failed; ",
                           "stability metrics are undefined (NA), not zero."), n_boot),
            call. = FALSE)
  } else if (failed_refits > 0L) {
    warning(sprintf(paste0("assess_program_stability(): %d of %d subsample refits failed; ",
                           "conditional_recovery_freq is computed over the %d successful ",
                           "refits (see replicate_log)."),
                    failed_refits, n_boot, valid_refits), call. = FALSE)
  }
  stab <- do.call(rbind, lapply(prog_ids, function(p) {
    r <- rec[[p]]
    data.frame(
      program_id = p,
      cell_type = prog_tab$cell_type[match(p, prog_ids)],
      n_hit = r$hit,
      n_eligible = r$eligible,
      recovery_freq = if (valid_refits > 0L) r$hit / n_boot else NA_real_,
      conditional_recovery_freq = if (valid_refits > 0L) r$hit / valid_refits else NA_real_,
      median_score_cor = if (length(r$score_cor)) stats::median(r$score_cor) else NA_real_,
      median_loading_cos = if (length(r$loading_cos)) stats::median(r$loading_cos) else NA_real_,
      sign_consistency = if (r$hit > 0L) r$sign_consistent / r$hit else NA_real_,
      top_gene_jaccard = if (length(r$top_jaccard)) stats::median(r$top_jaccard) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  repl <- do.call(rbind, repl_log)
  fit$stability <- list(per_program = stab, n_boot = n_boot,
                        valid_refits = valid_refits, failed_refits = failed_refits,
                        replicate_log = repl,
                        match_method = match_method, sim_threshold = sim_threshold)
  fit
}
