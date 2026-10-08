# cellprograms: first-stage EBMF backend (flashier) + canonicalization
#
# Implements plan Sections 9-12:
#   Y_c = Z_c W_c' + E_c per cell type, flashier backend,
#   RMS-scale + sign canonicalization that preserves reconstruction exactly.

#' Fit cell-type-specific programs with flashier
#'
#' @param x A `cell_program_data` object.
#' @param loading_prior Gene-side EBNM prior family: "point_normal" (default),
#'   "point_laplace", or "unimodal".
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
                                  covariates = NULL, seed = 1,
                                  flash_control = list()) {
  stopifnot(inherits(x, "cell_program_data"))
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
    if (!is.null(covariates)) Y <- .residualize(Y, x$sample_metadata, covariates, ct)
    if (center) Y <- scale(Y, center = TRUE, scale = scale)
    else if (scale) Y <- scale(Y, center = FALSE, scale = TRUE)

    # Drop genes that are identically zero after preprocessing (silent in this
    # cell type, or constant and centered away): flashier rejects all-zero
    # columns and they carry no factor signal.
    keep_cols <- colSums(abs(Y)) > 0
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
    fit <- tryCatch(
      do.call(flashier::flash, args),
      error = function(e) {
        # point_laplace's nlm solver can fail numerically on some pseudobulk
        # profiles; degrade gracefully to point_normal and record it.
        if (identical(prior_fn, ebnm::ebnm_point_laplace)) {
          cli::cli_warn("EBNM solver failed with point_laplace in cell type {.val {ct}}; retrying with point_normal.")
          args$ebnm_fn <- list(ebnm::ebnm_normal, ebnm::ebnm_point_normal)
          do.call(flashier::flash, args)
        } else {
          cli::cli_abort("flashier failed in cell type {.val {ct}}: {conditionMessage(e)}")
        }
      }
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    fits[[ct]] <- list(
      flash = fit,
      Y_used = Y,
      genes = genes,
      elapsed_secs = elapsed
    )
  }

  x$fits <- fits
  x$preprocessing <- list(
    centered = center, scaled = scale,
    feature_sets = feature_sets,
    loading_prior = loading_prior,
    max_factors = max_factors, var_type = var_type,
    backfit = backfit, nullcheck = nullcheck,
    covariates = covariates
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

# Residualize Y (obs x genes) on sample-level covariates, intercept included.
# lm.fit pivots, so covariate levels absent from this cell type's donors are
# dropped rather than erroring; only too few donors or NA covariates abort.
.residualize <- function(Y, meta, covariates, ct) {
  f <- if (inherits(covariates, "formula")) covariates
       else stats::reformulate(covariates)
  vars <- all.vars(f)
  absent <- setdiff(vars, names(meta))
  if (length(absent)) {
    cli::cli_abort("{.arg covariates} not in {.field sample_metadata}: {.val {absent}}.")
  }
  m <- meta[match(rownames(Y), meta$observation_id), vars, drop = FALSE]
  if (anyNA(m)) {
    cli::cli_abort("Missing covariate values among donors of cell type {.val {ct}}; impute or drop them first.")
  }
  X <- stats::model.matrix(f, m)
  if (nrow(X) <= qr(X)$rank + 1L) {
    cli::cli_abort("Too few donors ({nrow(X)}) in {.val {ct}} for {qr(X)$rank} covariate column(s).")
  }
  R <- stats::lm.fit(X, Y)$residuals
  dimnames(R) <- dimnames(Y)
  R
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
