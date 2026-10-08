## fit-ebmf.R — stage-1 per-cell-type EBMF via flashier,
## with Level-1 (shared priors) and Level-2 (initialization borrowing) hooks.

## ---------------------------------------------------------------------------
## EBNM function builders
## ---------------------------------------------------------------------------

.ebnm_family <- function(prior) {
  switch(
    prior,
    point_normal  = ebnm::ebnm_point_normal,
    point_laplace = ebnm::ebnm_point_laplace,
    normal        = ebnm::ebnm_normal,
    .stopf("Unknown prior family '%s'.", prior)
  )
}

## Closure with a fixed prior g (used by L1 shared-prior refits).
## flashier's validator calls ebnm_fn(x, s, g_init = NULL, fix_g = FALSE, output);
## we deliberately ignore those and always use the fixed pooled prior.
.make_fixed_prior_fn <- function(g) {
  force(g)
  function(x, s, g_init, fix_g, output) {
    ebnm::ebnm(x, s, prior_family = "normal_scale_mixture",
               g_init = g, fix_g = TRUE, output = output)
  }
}

## Pool fitted per-factor priors by mixture averaging.
## Each element of g_list is a normalmix (the gene-side prior of one factor);
## the pooled prior is the equally weighted average of these mixtures.
.pool_priors <- function(g_list) {
  g_list <- Filter(Negate(is.null), g_list)
  if (length(g_list) == 0L) return(NULL)
  w   <- unlist(lapply(g_list, function(g) g$pi))
  mns <- unlist(lapply(g_list, function(g) g$mean))
  sds <- unlist(lapply(g_list, function(g) g$sd))
  w <- w / length(g_list)
  keep <- w > 1e-12
  ashr::normalmix(pi = w[keep], mean = mns[keep], sd = sds[keep])
}

## Summary moments of a normalmix prior, for convergence checks and reporting.
.prior_moments <- function(g) {
  slab <- g$sd > 1e-8
  c(
    slab_weight = sum(g$pi[slab]),
    slab_sd     = if (any(slab)) sqrt(sum(g$pi[slab] * (g$sd[slab]^2 + g$mean[slab]^2)) /
                                        max(sum(g$pi[slab]), 1e-12)) else 0
  )
}

## Extract the fitted gene-side (mode 2) prior of each retained factor.
.factor_priors <- function(f) {
  gs <- f$flash_fit$g
  if (is.null(gs)) return(list())
  lapply(gs, function(g) if (length(g) >= 2L) g[[2L]] else NULL)
}

## ---------------------------------------------------------------------------
## Preprocessing
## ---------------------------------------------------------------------------

.select_features <- function(x, n_variable = NULL) {
  if (is.list(x$features)) return(x$features)
  out <- vector("list", length(x$cell_types))
  names(out) <- x$cell_types
  for (ct in x$cell_types) {
    m <- x$matrices[[ct]]
    gv <- apply(m, 2L, stats::var, na.rm = TRUE)
    gv[is.na(gv)] <- 0
    keep <- colnames(m)[gv > .Machine$double.eps]
    if (identical(x$features, "variable") && !is.null(n_variable) &&
        length(keep) > n_variable) {
      keep <- colnames(m)[order(-gv)][seq_len(n_variable)]
    }
    out[[ct]] <- keep
  }
  out
}

.preprocess <- function(x, feature_sets, center, scale) {
  mats <- vector("list", length(x$cell_types)); names(mats) <- x$cell_types
  centers <- vector("list", length(x$cell_types)); names(centers) <- x$cell_types
  scales  <- vector("list", length(x$cell_types)); names(scales)  <- x$cell_types
  for (ct in x$cell_types) {
    m <- x$matrices[[ct]][, feature_sets[[ct]], drop = FALSE]
    cm <- colMeans(m, na.rm = TRUE)
    if (center) m <- sweep(m, 2L, cm, `-`)
    cs <- apply(m, 2L, stats::sd, na.rm = TRUE)
    cs[cs < .Machine$double.eps] <- 1
    if (scale) m <- sweep(m, 2L, cs, `/`)
    mats[[ct]] <- m
    centers[[ct]] <- if (center) cm else stats::setNames(rep(0, length(cm)), names(cm))
    scales[[ct]]  <- if (scale) cs else stats::setNames(rep(1, length(cs)), names(cs))
  }
  list(matrices = mats, centers = centers, scales = scales)
}

## ---------------------------------------------------------------------------
## Single-cell-type fits
## ---------------------------------------------------------------------------

.fit_one_flash <- function(Y, score_prior, loading_prior, max_factors,
                           var_type, backfit, nullcheck, fixed_g = NULL,
                           verbose = 0L) {
  score_fn <- .ebnm_family(score_prior)
  load_fn  <- if (!is.null(fixed_g)) .make_fixed_prior_fn(fixed_g) else .ebnm_family(loading_prior)
  flashier::flash(
    Y,
    ebnm_fn     = list(score_fn, load_fn),
    var_type    = var_type,
    greedy_Kmax = max_factors,
    backfit     = backfit,
    nullcheck   = nullcheck,
    verbose     = verbose
  )
}

## L2: build a flashier initialization from donor loadings.
## donor_loadings must have gene row names; genes absent from the donor
## contribute zero (union representation, zero-padded).
.build_init <- function(Y_target, donor_loadings) {
  genes <- colnames(Y_target)
  F_init <- matrix(0, nrow = length(genes), ncol = ncol(donor_loadings),
                   dimnames = list(genes, colnames(donor_loadings)))
  common <- intersect(genes, rownames(donor_loadings))
  if (length(common) == 0L) {
    .stopf("init_from: no overlapping genes between donor loadings and target matrix.")
  }
  F_init[common, ] <- donor_loadings[common, , drop = FALSE]
  ## initialize scores by ridge regression of Y on donor loadings
  FtF <- crossprod(F_init) + 1e-6 * diag(ncol(F_init))
  L_init <- Y_target %*% F_init %*% MASS_ginv(FtF)
  list(L_init, F_init)
}

## Minimal Moore-Penrose inverse without importing MASS at namespace level
MASS_ginv <- function(X) {
  s <- svd(X)
  d <- s$d
  d[d > 1e-10] <- 1 / d[d > 1e-10]
  s$v %*% diag(d, nrow = length(d)) %*% t(s$u)
}

.fit_one_with_init <- function(Y, init, score_prior, loading_prior,
                               extra_Kmax, var_type, nullcheck, fix_borrowed,
                               verbose = 0L) {
  ebnm_list <- list(.ebnm_family(score_prior), .ebnm_family(loading_prior))
  fl <- flashier::flash_init(Y, var_type = var_type) |>
    flashier::flash_factors_init(init, ebnm_fn = ebnm_list)
  if (fix_borrowed) {
    fl <- flashier::flash_factors_fix(fl, kset = seq_len(ncol(init[[2L]])),
                                      which_dim = "F")
  }
  if (extra_Kmax > 0L) {
    fl <- fl |> flashier::flash_greedy(Kmax = extra_Kmax, ebnm_fn = ebnm_list,
                                       verbose = 0L)
  }
  fl <- fl |> flashier::flash_backfit(verbose = 0L)
  if (fix_borrowed) {
    fl <- flashier::flash_factors_unfix(fl, kset = seq_len(ncol(init[[2L]])),
                                        which_dim = "F")
  }
  if (nullcheck) {
    fl <- fl |> flashier::flash_nullcheck(verbose = 0L)
  }
  fl
}

## ---------------------------------------------------------------------------
## Extraction
## ---------------------------------------------------------------------------

.extract_fit <- function(f, ct, Y = NULL) {
  K <- f$n_factors
  ids <- .program_ids(ct, K)
  Z <- f$L_pm; W <- f$F_pm
  if (K == 0L || is.null(Z) || is.null(W) || ncol(Z) == 0L) {
    ## 0-factor fit: return properly dimensioned empty matrices so downstream
    ## extractors (ncol, rownames, program_scores) never see NULL
    n <- if (!is.null(Y)) nrow(Y) else 0L
    g <- if (!is.null(Y)) ncol(Y) else 0L
    Z <- matrix(numeric(0), n, 0L, dimnames = list(rownames(Y), character(0)))
    W <- matrix(numeric(0), g, 0L, dimnames = list(colnames(Y), character(0)))
    return(list(scores = Z, loadings = W, scores_sd = Z, loadings_sd = W,
                lfsr = NULL, pve = numeric(0), residuals_sd = f$residuals_sd,
                elbo = f$elbo, flash = f))
  }
  colnames(Z) <- ids; rownames(W) <- rownames(f$F_pm) %||% rownames(W)
  colnames(W) <- ids
  list(
    scores      = Z,
    loadings    = W,
    scores_sd   = f$L_psd,
    loadings_sd = f$F_psd,
    lfsr        = f$F_lfsr,
    pve         = f$pve,
    residuals_sd = f$residuals_sd,
    elbo        = f$elbo,
    flash       = f
  )
}

## ---------------------------------------------------------------------------
## Main entry point
## ---------------------------------------------------------------------------

#' Fit per-cell-type empirical-Bayes matrix factorization
#'
#' Stage-1 local factorization. Each cell type is fit independently;
#' expression matrices are never concatenated before factorization.
#'
#' @param x A `cell_program_data` object.
#' @param loading_prior Gene-side EBNM prior: `"point_normal"` (default),
#'   `"point_laplace"`, or `"normal"`.
#' @param score_prior Sample-side EBNM prior (default `"normal"`).
#' @param max_factors Maximum candidate factors per cell type (greedy Kmax).
#' @param center Gene-center each cell-type matrix (default TRUE).
#' @param scale Unit-variance scale genes (default FALSE).
#' @param var_type flashier variance structure; 2 = per-gene (default).
#' @param backfit Backfit after greedy fitting (default TRUE).
#' @param nullcheck Null-check factors after fitting (default TRUE).
#' @param n_variable Optional cap on per-cell-type features (top by variance)
#'   when `features = "variable"`.
#' @param share_prior Level-1 borrowing: `"none"` (default), `"global"`, or
#'   `"groups"` (pool EBNM loading priors across cell types / within
#'   `prior_groups`, then refit with the pooled prior fixed).
#' @param prior_groups Named list mapping group labels to cell-type vectors;
#'   used when `share_prior = "groups"`.
#' @param prior_max_iter Maximum outer empirical-Bayes iterations.
#' @param prior_tol Convergence tolerance on pooled prior moments.
#' @param init_from Level-2 borrowing: named list, one entry per target cell
#'   type, each a list with elements `fit` (a donor `cell_program_fit`),
#'   `cell_type` (donor cell type), and optionally `programs` (program IDs;
#'   default all). Donor loadings initialize the target fit and may shrink
#'   to zero during backfitting.
#' @param fix_borrowed If TRUE, borrowed loadings are fixed during backfit
#'   (reference-projection mode). Default FALSE.
#' @param seed Random seed (flashier greedy initialization is randomized).
#' @param verbose Verbosity level passed to flashier.
#' @return A `cell_program_fit` object.
#' @export
fit_celltype_programs <- function(x,
                                  loading_prior = c("point_normal", "point_laplace", "normal"),
                                  score_prior = "normal",
                                  max_factors = 30L,
                                  center = TRUE,
                                  scale = FALSE,
                                  var_type = 2L,
                                  backfit = TRUE,
                                  nullcheck = TRUE,
                                  n_variable = NULL,
                                  share_prior = c("none", "global", "groups"),
                                  prior_groups = NULL,
                                  prior_max_iter = 5L,
                                  prior_tol = 1e-3,
                                  init_from = NULL,
                                  fix_borrowed = FALSE,
                                  seed = NULL,
                                  verbose = 0L) {
  validate_cell_program_data(x)
  loading_prior <- match.arg(loading_prior)
  share_prior <- match.arg(share_prior)
  if (!is.null(seed)) set.seed(seed)

  feature_sets <- .select_features(x, n_variable)
  pre <- .preprocess(x, feature_sets, center, scale)
  mats <- pre$matrices

  ## ---- Level 2: build initializations ------------------------------------
  inits <- vector("list", length(x$cell_types)); names(inits) <- x$cell_types
  l2_report <- list()
  if (!is.null(init_from)) {
    for (ct in names(init_from)) {
      if (!ct %in% x$cell_types) .stopf("init_from: unknown target cell type '%s'.", ct)
      spec <- init_from[[ct]]
      donor_fit <- spec$fit
      donor_ct  <- spec$cell_type
      if (!inherits(donor_fit, "cell_program_fit")) {
        .stopf("init_from[['%s']]$fit must be a cell_program_fit.", ct)
      }
      donor_W <- donor_fit$loadings[[donor_ct]]
      if (!is.null(spec$programs)) donor_W <- donor_W[, spec$programs, drop = FALSE]
      inits[[ct]] <- .build_init(mats[[ct]], donor_W)
      l2_report[[ct]] <- list(
        donor_cell_type = donor_ct,
        n_borrowed = ncol(donor_W),
        n_overlap_genes = sum(rownames(donor_W) %in% colnames(mats[[ct]]))
      )
    }
  }

  ## ---- Stage-1 fits --------------------------------------------------------
  fit_all <- function(fixed_priors = NULL) {
    fits <- vector("list", length(x$cell_types)); names(fits) <- x$cell_types
    for (ct in x$cell_types) {
      fg <- if (!is.null(fixed_priors)) fixed_priors[[ct]] else NULL
      if (!is.null(inits[[ct]])) {
        extra <- max(0L, max_factors - ncol(inits[[ct]][[2L]]))
        fits[[ct]] <- .fit_one_with_init(mats[[ct]], inits[[ct]], score_prior,
                                         loading_prior, extra, var_type,
                                         nullcheck, fix_borrowed, verbose)
      } else {
        fits[[ct]] <- .fit_one_flash(mats[[ct]], score_prior, loading_prior,
                                     max_factors, var_type, backfit, nullcheck,
                                     fixed_g = fg, verbose = verbose)
      }
    }
    fits
  }

  fits <- fit_all()

  ## ---- Level 1: shared-prior outer EB loop --------------------------------
  l1_report <- NULL
  if (share_prior != "none") {
    groups <- if (share_prior == "global") {
      list(all = x$cell_types)
    } else {
      if (is.null(prior_groups)) .stopf("`prior_groups` required when share_prior = 'groups'.")
      prior_groups
    }
    trace <- list()
    prev_moments <- NULL
    for (iter in seq_len(prior_max_iter)) {
      pooled <- vector("list", length(x$cell_types)); names(pooled) <- x$cell_types
      moments <- c()
      for (grp in groups) {
        gs <- unlist(lapply(fits[grp], .factor_priors), recursive = FALSE)
        g_pool <- .pool_priors(gs)
        if (is.null(g_pool)) next
        for (ct in grp) pooled[[ct]] <- g_pool
        moments <- rbind(moments, .prior_moments(g_pool))
      }
      trace[[iter]] <- moments
      if (!is.null(prev_moments) &&
          max(abs(moments - prev_moments)) < prior_tol) break
      prev_moments <- moments
      fits <- fit_all(fixed_priors = pooled)
    }
    l1_report <- list(
      mode = share_prior, groups = groups, n_iter = iter,
      pooled_priors = pooled, trace = trace,
      moments = stats::setNames(as.data.frame(do.call(rbind, trace)), NULL)
    )
  }

  ## ---- Extract -------------------------------------------------------------
  scores <- loadings <- list()
  scores_sd <- loadings_sd <- lfsr <- list()
  for (ct in x$cell_types) {
    ex <- .extract_fit(fits[[ct]], ct, mats[[ct]])
    scores[[ct]] <- ex$scores;       loadings[[ct]] <- ex$loadings
    scores_sd[[ct]] <- ex$scores_sd; loadings_sd[[ct]] <- ex$loadings_sd
    lfsr[[ct]] <- ex$lfsr
  }

  programs <- .build_program_table(fits, scores, loadings, lfsr, x$cell_types)

  structure(
    list(
      input = list(
        matrices = mats,
        observation_ids = lapply(mats, rownames),
        sample_metadata = x$sample_metadata,
        feature_sets = feature_sets
      ),
      preprocessing = list(center = center, scale = scale,
                           centers = pre$centers, scales = pre$scales),
      fits = fits,
      programs = programs,
      scores = scores,
      loadings = loadings,
      posteriors = list(scores_sd = scores_sd, loadings_sd = loadings_sd, lfsr = lfsr),
      canonicalization = list(applied = FALSE),
      stability = NULL,
      borrowing = list(level1 = l1_report,
                       level2 = if (length(l2_report)) list(init_map = l2_report) else NULL),
      alignment = NULL,
      integration = NULL,
      provenance = list(
        package_version = as.character(utils::packageVersion("cellprograms") %||% "0.1.0"),
        flashier_version = as.character(utils::packageVersion("flashier")),
        ebnm_version = as.character(utils::packageVersion("ebnm")),
        seed = seed, timestamp = Sys.time(),
        parameters = list(loading_prior = loading_prior, score_prior = score_prior,
                          max_factors = max_factors, center = center, scale = scale,
                          var_type = var_type, backfit = backfit, nullcheck = nullcheck,
                          share_prior = share_prior, fix_borrowed = fix_borrowed)
      )
    ),
    class = "cell_program_fit"
  )
}

.build_program_table <- function(fits, scores, loadings, lfsr, cell_types) {
  rows <- list()
  for (ct in cell_types) {
    K <- fits[[ct]]$n_factors
    if (K == 0L) next
    W <- loadings[[ct]]
    for (k in seq_len(K)) {
      w <- W[, k]
      rows[[length(rows) + 1L]] <- data.frame(
        program_id = colnames(W)[k],
        cell_type = ct,
        factor_index = k,
        pve = fits[[ct]]$pve[k] %||% NA_real_,
        sparsity_lfsr05 = if (!is.null(lfsr[[ct]])) mean(lfsr[[ct]][, k] < 0.05) else NA_real_,
        max_abs_loading = max(abs(w)),
        rms_loading = .rms(w),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(program_id = character(0), cell_type = character(0)))
  }
  do.call(rbind, rows)
}

#' @export
print.cell_program_fit <- function(x, ...) {
  cat("cell_program_fit\n")
  for (ct in names(x$fits)) {
    cat(sprintf("    %-12s %2d programs  (pve: %s)\n",
                ct, x$fits[[ct]]$n_factors,
                paste(round(x$fits[[ct]]$pve, 2), collapse = ", ")))
  }
  if (!is.null(x$borrowing$level1)) {
    cat("  L1 shared prior:", x$borrowing$level1$mode,
        sprintf("(%d EB iterations)\n", x$borrowing$level1$n_iter))
  }
  if (!is.null(x$borrowing$level2)) {
    cat("  L2 init borrowing for:", paste(names(x$borrowing$level2$init_map), collapse = ", "), "\n")
  }
  if (!is.null(x$integration)) cat("  integration:", x$integration$method, "\n")
  invisible(x)
}

#' @export
summary.cell_program_fit <- function(object, ...) {
  object$programs
}
