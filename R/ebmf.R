# cellprograms: first-stage EBMF backend (flashier) + canonicalization
#
# Implements plan Sections 9-12:
#   Y_c = Z_c W_c' + E_c per cell type, flashier backend,
#   RMS-scale + sign canonicalization that preserves reconstruction exactly.

#' Fit cell-type-specific programs with flashier
#'
#' @param x A `cell_program_data` object.
#' @param loading_prior Gene-side EBNM prior family: "point_laplace" (default),
#'   "point_normal", or "unimodal".
#' @param max_factors Maximum candidate factors per cell type (greedy_Kmax).
#'   Default 10: benchmarked on COMBAT pseudobulk (32-config factorial),
#'   K=10 dominated all prior/variance/backfit combinations; K >= 30
#'   overfits and collapses disease-signal recovery.
#' @param center Center genes within each cell type (default TRUE).
#' @param scale Unit-variance scale genes (default FALSE; not recommended).
#' @param features Feature selection mode: "all", "variable", or a named list
#'   of gene ID vectors per cell type.
#' @param n_variable_genes Number of most-variable genes when
#'   `features = "variable"`.
#' @param var_type flashier residual-variance type: 0 = constant, 1 = one
#'   variance per row/donor (default), 2 = per-gene (unstable on
#'   log-normalized pseudobulk), or `c(1, 2)` = Kronecker rank-one structure
#'   s_ij = a_i * b_j (most flexible; substantially slower).
#' @param backfit Run flashier backfitting (default TRUE).
#' @param nullcheck Run flashier nullcheck (default TRUE).
#' @param covariates Optional sample-level covariates to regress out of each
#'   cell type's pseudobulk (genes x samples fit separately per cell type, over
#'   the donors observed in that cell type) before factorization: a character
#'   vector of `sample_metadata` column names or a one-sided formula, e.g.
#'   `~ age + sex + batch`. Programs are then fit to the residuals. Do not
#'   include the outcome you plan to test. Variable-gene selection is done
#'   before adjustment.
#' @param covariate_mode How `covariates` are handled: `"residualize"`
#'   (default; hard adjustment: regress out of `Y` before factorization) or
#'   `"fixed"` (soft, SOFA-style: the covariate design columns enter flashier as
#'   fixed sample-side factors with unshrunk gene effects, and the K free
#'   programs are fit jointly with them; free programs can still share
#'   variance with a covariate. Requires `center = TRUE`, ignores
#'   `flash_control`). Fixed columns are not returned as programs;
#'   `Y_used` is `Y` minus the fitted covariate part.
#' @param seed Random seed for reproducibility.
#' @param flash_control Optional list overriding flashier::flash arguments.
#'
#' @return The input object augmented with class `cell_program_fit` and
#'   components `fits`, `programs`, `scores`, `loadings`, `preprocessing`,
#'   `provenance`.
#'
#' @export
fit_celltype_programs <- function(x, loading_prior = "point_laplace",
                                  max_factors = 10, center = TRUE,
                                  scale = FALSE, features = "all",
                                  n_variable_genes = 2000, var_type = 1,
                                  backfit = TRUE, nullcheck = TRUE,
                                  covariates = NULL,
                                  covariate_mode = c("residualize", "fixed"),
                                  seed = 1, flash_control = list()) {
  stopifnot(inherits(x, "cell_program_data"))
  covariate_mode <- match.arg(covariate_mode)
  if (covariate_mode == "fixed" && !is.null(covariates) && !center) {
    cli::cli_abort("{.code covariate_mode = \"fixed\"} requires {.code center = TRUE}.")
  }
  if (!requireNamespace("flashier", quietly = TRUE) ||
      !requireNamespace("ebnm", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg flashier} and {.pkg ebnm} are required.")
  }
  if (!is.numeric(var_type) || length(var_type) > 2L ||
      !all(var_type %in% c(0, 1, 2)) ||
      (length(var_type) == 2L && !setequal(var_type, c(1, 2)))) {
    cli::cli_abort("{.arg var_type} must be 0, 1, 2, or c(1, 2) (Kronecker).")
  }

  prior_fn <- switch(
    loading_prior,
    point_normal = ebnm::ebnm_point_normal,
    point_laplace = ebnm::ebnm_point_laplace,
    unimodal = ebnm::ebnm_unimodal,
    cli::cli_abort("Unknown {.arg loading_prior}: {.val {loading_prior}}.")
  )

  feature_sets <- .select_features(x, features, n_variable_genes)

  fits <- list()
  for (ct in x$cell_types) {
    Y <- as.matrix(x$pseudobulk[[ct]])
    genes <- feature_sets[[ct]]
    Y <- Y[, genes, drop = FALSE]
    if (!is.null(covariates) && covariate_mode == "residualize") {
      Y <- .residualize(Y, x$sample_metadata, covariates, ct)
    }
    if (center) Y <- scale(Y, center = TRUE, scale = scale)
    else if (scale) Y <- scale(Y, center = FALSE, scale = TRUE)
    Xc <- if (!is.null(covariates) && covariate_mode == "fixed") {
      .cov_design(rownames(Y), x$sample_metadata, covariates, ct, intercept = FALSE)
    }

    # Drop genes that are identically zero after preprocessing (silent in this
    # cell type, or constant and centered away): flashier rejects all-zero
    # columns and they carry no factor signal.
    cs <- colSums(abs(Y))
    keep_cols <- cs > 1e-10 * max(cs, 1)   # tolerance: centered constant genes leave ~1e-16 dust
    if (!any(keep_cols)) {
      # No informative genes in this cell type (e.g. a pseudobulk profile that
      # is constant across samples): record an empty fit (K=0) instead of
      # failing the whole analysis. Downstream, canonicalize_programs() and
      # the scores/loadings extraction handle K=0 cell types.
      cli::cli_inform("All genes are degenerate (zero after preprocessing) in cell type {.val {ct}}; recording K=0 fit.")
      fits[[ct]] <- list(
        flash = list(L_pm = NULL, F_pm = NULL),
        Y_used = Y,
        genes = character(0),
        elapsed_secs = 0,
        skipped = TRUE
      )
      next
    }
    if (!all(keep_cols)) {
      cli::cli_inform("Dropping {sum(!keep_cols)} zero-variance gene(s) in cell type {.val {ct}}.")
      Y <- Y[, keep_cols, drop = FALSE]
      genes <- genes[keep_cols]
    }

    args <- modifyList(list(
      data = Y,
      ebnm_fn = list(ebnm::ebnm_normal, prior_fn),
      var_type = var_type,
      greedy_Kmax = max_factors,
      backfit = backfit,
      nullcheck = nullcheck,
      verbose = 0L
    ), flash_control)

    t0 <- Sys.time()
    set.seed(seed)
    run_flash <- function(args) {
      if (is.null(Xc)) do.call(flashier::flash, args) else .flash_fixed(args, Xc)
    }
    # Fallback ladder for numerical failures in flashier: point_laplace's nlm
    # solver can fail on some pseudobulk profiles (-> point_normal), and the
    # nullcheck can hit a NaN ELBO when K approaches the residual degrees of
    # freedom (-> rerun without nullcheck). Every fallback is warned about and
    # recorded in fits[[ct]]$fallback.
    attempts <- list(list(label = "none", args = args))
    if (identical(prior_fn, ebnm::ebnm_point_laplace)) {
      a2 <- args; a2$ebnm_fn <- list(ebnm::ebnm_normal, ebnm::ebnm_point_normal)
      attempts[[length(attempts) + 1L]] <- list(label = "point_normal", args = a2)
    }
    if (isTRUE(args$nullcheck)) {
      for (at in attempts) {
        a3 <- at$args; a3$nullcheck <- FALSE
        attempts[[length(attempts) + 1L]] <- list(
          label = if (at$label == "none") "no_nullcheck" else paste0(at$label, "+no_nullcheck"),
          args = a3)
      }
    }
    fit <- NULL; fallback <- "none"; first_err <- NULL
    for (at in attempts) {
      set.seed(seed)                       # every attempt starts from the same RNG state
      fit <- tryCatch(run_flash(at$args), error = function(e) {
        if (is.null(first_err)) first_err <<- e
        NULL
      })
      if (!is.null(fit) && !.fit_ok(fit)) {   # silent NaN ELBO is a failure too
        if (is.null(first_err)) first_err <<- simpleError("flashier returned a non-finite ELBO")
        fit <- NULL
      }
      if (!is.null(fit)) { fallback <- at$label; break }
    }
    if (is.null(fit)) {
      cli::cli_abort("flashier failed in cell type {.val {ct}}: {conditionMessage(first_err)}")
    }
    if (fallback != "none") {
      cli::cli_warn("flashier needed fallback {.val {fallback}} in cell type {.val {ct}}.")
    }
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!is.null(fit$fixed_L)) Y <- Y - fit$fixed_L %*% t(fit$fixed_F)
    fits[[ct]] <- list(
      flash = fit,
      Y_used = Y,
      genes = genes,
      elapsed_secs = elapsed,
      fallback = fallback
    )
  }

  x$fits <- fits
  x$preprocessing <- list(
    centered = center, scaled = scale,
    feature_sets = feature_sets,
    loading_prior = loading_prior,
    max_factors = max_factors, var_type = var_type,
    backfit = backfit, nullcheck = nullcheck,
    covariates = covariates, covariate_mode = covariate_mode
  )
  x$provenance <- list(
    seed = seed,
    timestamp = Sys.time(),
    r_version = R.version.string,
    flashier_version = as.character(utils::packageVersion("flashier")),
    ebnm_version = as.character(utils::packageVersion("ebnm"))
  )
  x$canonicalization <- NULL
  x$stability <- NULL

  class(x) <- c("cell_program_fit", "cell_program_data", "list")
  x
}

# A flashier fit is usable if its ELBO (when reported) is finite.
.fit_ok <- function(fit) is.null(fit$elbo) || all(is.finite(fit$elbo))

# Design matrix of sample-level covariates for the donors in `ids`.
# Terms constant across this cell type's donors (e.g. a site with no donors
# here) cannot be adjusted for and model.matrix() errors on one-level factors,
# so they are dropped for this cell type only. NULL if nothing is left.
# intercept = FALSE: centered, unit-variance, full-rank columns (fixed mode).
.cov_design <- function(ids, meta, covariates, ct, intercept = TRUE) {
  f <- if (inherits(covariates, "formula")) covariates
       else stats::reformulate(covariates)
  vars <- all.vars(f)
  absent <- setdiff(vars, names(meta))
  if (length(absent)) {
    cli::cli_abort("{.arg covariates} not in {.field sample_metadata}: {.val {absent}}.")
  }
  m <- meta[match(ids, meta$observation_id), vars, drop = FALSE]
  if (anyNA(m)) {
    cli::cli_abort("Missing covariate values among donors of cell type {.val {ct}}; impute or drop them first.")
  }
  m[] <- lapply(m, function(v) if (is.factor(v)) droplevels(v) else v)   # unused levels -> all-zero dummies
  const <- vapply(vars, function(v) length(unique(m[[v]])) < 2L, logical(1L))
  if (all(const)) {
    cli::cli_inform("Covariate(s) {.val {vars}} constant in cell type {.val {ct}}; not adjusted.")
    return(NULL)
  }
  if (any(const)) {
    cli::cli_inform("Covariate(s) {.val {vars[const]}} constant in cell type {.val {ct}}; dropped for this cell type.")
    f <- stats::drop.terms(stats::terms(f), which(attr(stats::terms(f), "term.labels") %in% vars[const]),
                           keep.response = FALSE)
    m <- m[, vars[!const], drop = FALSE]
  }
  X <- stats::model.matrix(f, m)
  rk <- qr(X)$rank
  if (nrow(X) - rk < 10L) {
    cli::cli_abort("Too few donors ({nrow(X)}) in {.val {ct}} for {rk} covariate column(s): need >= 10 residual degrees of freedom.")
  }
  hat <- rowSums(qr.Q(qr(X))[, seq_len(rk), drop = FALSE]^2)
  if (any(hat > 0.99)) {
    cli::cli_warn("{sum(hat > 0.99)} donor(s) in cell type {.val {ct}} have leverage ~1 under the covariate model (e.g. a level with a single donor): their residuals are ~0.")
  }
  if (!intercept) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    X <- scale(X)
    q <- qr(X)
    X <- X[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
    attr(X, "scaled:center") <- attr(X, "scaled:scale") <- NULL
  }
  X
}

# Residualize Y (obs x genes) on sample-level covariates, intercept included.
.residualize <- function(Y, meta, covariates, ct) {
  X <- .cov_design(rownames(Y), meta, covariates, ct)
  if (is.null(X)) return(Y)
  R <- stats::lm.fit(X, Y)$residuals
  dimnames(R) <- dimnames(Y)
  R
}

# Soft adjustment (SOFA-style): covariate columns enter flashier as FIXED
# sample-side factors with unshrunk (wide fixed-scale normal) gene effects,
# initialised at OLS; free programs are then added greedily and everything but
# the fixed sample-side columns is backfit. Returns a flash-like list holding
# only the free factors, plus fixed_L / fixed_F for the covariate part.
.flash_fixed <- function(args, Xc) {
  Y <- args$data; q <- ncol(Xc)
  F0 <- t(solve(crossprod(Xc), crossprod(Xc, Y)))
  fl <- flashier::flash_init(Y, var_type = args$var_type) |>
    flashier::flash_set_verbose(args$verbose) |>
    flashier::flash_factors_init(
      list(Xc, F0),
      ebnm_fn = flashier::flash_ebnm(prior_family = "normal", mode = 0, scale = 10)) |>
    flashier::flash_factors_fix(kset = seq_len(q), which_dim = "loadings") |>
    flashier::flash_greedy(Kmax = args$greedy_Kmax, ebnm_fn = args$ebnm_fn)
  free <- function(fl) setdiff(seq_len(fl$n_factors), seq_len(q))
  if (args$backfit) fl <- flashier::flash_backfit(fl, verbose = args$verbose)
  if (args$nullcheck && length(free(fl))) {
    fl <- flashier::flash_nullcheck(fl, kset = free(fl), verbose = args$verbose)
  }
  k <- free(fl)
  out <- list(
    n_factors = length(k), elbo = fl$elbo, residuals_sd = fl$residuals_sd,
    L_pm = fl$L_pm[, k, drop = FALSE], F_pm = fl$F_pm[, k, drop = FALSE],
    F_lfsr = if (!is.null(fl$F_lfsr)) fl$F_lfsr[, k, drop = FALSE],
    fixed_L = Xc, fixed_F = fl$F_pm[, seq_len(q), drop = FALSE]
  )
  if (length(k) == 0L) { out$L_pm <- NULL; out$F_pm <- NULL }
  out
}

# Feature selection: independent per cell type; never uses outcomes.
.select_features <- function(x, features, n_variable_genes) {
  if (is.list(features)) {
    missing_ct <- setdiff(x$cell_types, names(features))
    if (length(missing_ct)) {
      cli::cli_abort("{.arg features} list lacks cell type(s): {.val {missing_ct}}.")
    }
    return(lapply(x$cell_types, function(ct) {
      g <- features[[ct]]
      avail <- colnames(x$pseudobulk[[ct]])
      unknown <- setdiff(g, avail)
      if (length(unknown)) {
        cli::cli_abort("{length(unknown)} requested gene(s) absent from {.val {ct}}.")
      }
      g
    }) |> stats::setNames(x$cell_types))
  }
  out <- list()
  for (ct in x$cell_types) {
    genes <- colnames(x$pseudobulk[[ct]])
    if (identical(features, "all")) {
      out[[ct]] <- genes
    } else if (identical(features, "variable")) {
      v <- apply(as.matrix(x$pseudobulk[[ct]]), 2, stats::var)
      k <- min(n_variable_genes, length(v))
      out[[ct]] <- names(sort(v, decreasing = TRUE))[seq_len(k)]
    } else {
      cli::cli_abort("{.arg features} must be {.val all}, {.val variable}, or a named list.")
    }
  }
  out
}

#' Canonicalize program scale and sign
#'
#' RMS-scales sample scores to unit root-mean-square with the compensating
#' factor applied to gene loadings (reconstruction preserved exactly), and
#' orients each factor so its largest-|loading| gene is positive.
#'
#' @param fit A `cell_program_fit` object.
#'
#' @export
canonicalize_programs <- function(fit) {
  stopifnot(inherits(fit, "cell_program_fit"))
  canon <- list()
  for (ct in fit$cell_types) {
    f <- fit$fits[[ct]]$flash
    if (is.null(f$L_pm) || is.null(f$F_pm)) {
      # Skipped (fully degenerate) cell type: K=0.
      canon[[ct]] <- list(scale_multiplier = numeric(0),
                          sign_multiplier = numeric(0))
      next
    }
    Z <- as.matrix(f$L_pm)  # observations x K (sample-side posterior mean)
    W <- as.matrix(f$F_pm)  # genes x K (gene-side posterior mean)
    K <- ncol(Z)
    if (K == 0L) {
      canon[[ct]] <- list(scale_multiplier = numeric(0),
                          sign_multiplier = numeric(0))
      next
    }
    s_k <- sqrt(colMeans(Z^2))                       # RMS of each score column
    s_k[s_k == 0] <- 1                               # degenerate guard
    # sign: gene with largest |loading|; ties broken by gene name (deterministic)
    sign_k <- apply(W, 2, function(w) {
      am <- which(abs(w) == max(abs(w)))
      am <- am[order(fit$fits[[ct]]$genes[am], w[am])]
      sign(w[am[1]])
    })
    sign_k[sign_k == 0] <- 1
    Zc <- sweep(Z, 2, s_k, "/")
    Wc <- sweep(W, 2, s_k, "*")
    Zc <- sweep(Zc, 2, sign_k, "*")
    Wc <- sweep(Wc, 2, sign_k, "*")
    canon[[ct]] <- list(
      Z = Zc, W = Wc,
      scale_multiplier = unname(s_k),
      sign_multiplier = unname(sign_k)
    )
  }
  fit$canonicalization <- canon
  # Populate the documented `scores` / `loadings` components with the
  # canonicalized matrices, labelled by observation / gene / program IDs.
  fit$scores <- list()
  fit$loadings <- list()
  for (ct in fit$cell_types) {
    Y_used <- fit$fits[[ct]]$Y_used
    genes_ct <- fit$fits[[ct]]$genes
    K <- if (is.null(canon[[ct]]$Z)) 0L else ncol(canon[[ct]]$Z)
    # NB: paste0(ct, "_", seq_len(0)) returns ct_" " (length 1), not
    # character(0) — paste() treats zero-length args as "". Guard K=0.
    prog_names <- if (K > 0L) paste0(ct, "_", seq_len(K)) else character(0)
    Z <- matrix(0, nrow = nrow(Y_used), ncol = K,
                dimnames = list(rownames(Y_used), prog_names))
    W <- matrix(0, nrow = length(genes_ct), ncol = K,
                dimnames = list(genes_ct, prog_names))
    if (K > 0L) {
      Z[] <- canon[[ct]]$Z
      W[] <- canon[[ct]]$W
    }
    fit$scores[[ct]] <- Z
    fit$loadings[[ct]] <- W
  }
  fit
}
