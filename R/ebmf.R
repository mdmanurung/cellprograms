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
#'   Must leave at least one residual degree of freedom per cell type:
#'   the fit aborts if `max_factors > n - q - 1`, where n is the number of
#'   donors for the cell type and q the number of covariate columns adjusted
#'   for (safety guard; never binds at n >= 50 with `max_factors <= 10`).
#'   Default 10: benchmarked on COMBAT pseudobulk (32-config factorial),
#'   K=10 dominated all prior/variance/backfit combinations; K >= 30
#'   overfits and collapses disease-signal recovery. For coordinated-program
#'   detection (`find_coordinated_programs()`), use `max_factors = 30`: it
#'   raised CMP membership F1 in the cmp-bench-v1.2 held-out benchmark.
#' @param center Center genes within each cell type (default TRUE).
#' @param scale Unit-variance scale genes (default FALSE; not recommended).
#' @param features Feature selection mode: "all", "variable", or a named list
#'   of gene ID vectors per cell type.
#' @param n_variable_genes Number of most-variable genes when
#'   `features = "variable"`.
#' @param exclude_genes Optional character vector of gene IDs removed from
#'   every cell type's feature universe BEFORE feature selection, so excluded
#'   genes do not consume `n_variable_genes` slots (with
#'   `features = "variable"`) and never enter the fit. Exclusion changes the
#'   feature panel and is NOT the default; it is intended for artifact-gene
#'   classes such as V(D)J or common-germline-variation genes -- see
#'   `artifact_gene_sets` and `artifact_genes()`. With a named-list
#'   `features`, excluded genes are dropped from the requested vectors after
#'   validation.
#' @param exclude_artifacts Convenience preset: when TRUE, union
#'   `artifact_genes()` (V(D)J, genotype-marker, sex and canonical-marker
#'   sets) over the pseudobulk gene universe into `exclude_genes`, so the
#'   preset blacklist never enters the fit and does not consume
#'   `n_variable_genes` slots. Default FALSE: exclusion changes the feature
#'   panel, so it is opt-in; explicit `exclude_genes` entries are kept.
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
#' @param feature_scale Declared measurement units of the pseudobulk
#'   expression values, recorded in the fit manifest and required for
#'   `project_programs_to_cells()` (the projection oracle refuses queries
#'   whose declared units differ from the source fit's). "unknown" (default)
#'   keeps the fit unprojectable; pass e.g. "logcounts" for
#'   log-normalized pseudobulk. Only used on the flashier.utils backend.
#' @param preprocessing_spec Optional frozen preprocessing specification from
#'   `freeze_preprocessing()` (MR-0). When supplied, the spec governs feature
#'   selection and adjustment: `features`, `n_variable_genes`,
#'   `exclude_genes`, `exclude_artifacts`, `covariates`, `covariate_mode`,
#'   `center`, and `scale` must all be left at their defaults (any explicit
#'   conflict is an error, so a refit can never silently drift from the
#'   frozen panel), and the spec's digest is recorded in the fit's
#'   preprocessing provenance. The spec's cell types must cover `x` and all
#'   frozen genes must be present in the dataset.
#' @param backend EBMF implementation: `"flashier"` (default; the
#'   flashier.utils R backend, with a nullcheck fallback ladder) or
#'   `"cebmf_torch"` (the covariate-EBMF torch backend, which fails loudly
#'   instead of falling back). Both backends share preprocessing, program
#'   extraction, and the output contract. The torch backend runs in a pinned
#'   Python environment: set options(cellprograms.cebmf_python_env = ...) to
#'   point at an existing environment, otherwise a basilisk environment is
#'   used.
#' @param var_floor Residual SD handling: `"none"` (default; residual SDs are
#'   estimated as usual), `"vignette"` (fixed at the flashier single-cell
#'   vignette's S = sd(log(Poisson(1/n) + 1)), drawn once per cell type with a
#'   fixed local seed, so the value is deterministic and independent of
#'   `seed`), or a positive numeric scalar S (residual SDs fixed at S). Fixed
#'   S binds on the `cebmf_torch` backend. On the `flashier` backend the
#'   value is passed through but does not currently bind (CAND-13: fits with
#'   and without a fixed scalar S are bit-identical at the pinned flashier
#'   version); a fail-closed remediation is pending.
#'
#' @return The input object augmented with class `cell_program_fit` and
#'   components `fits`, `programs`, `scores`, `loadings`, `preprocessing`,
#'   `provenance`.
#'
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable",
#'                              n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' program_summary(fit)
#' @export
fit_celltype_programs <- function(x, loading_prior = "point_laplace",
                                  max_factors = 10, center = TRUE,
                                  scale = FALSE, features = "all",
                                  n_variable_genes = 2000, exclude_genes = NULL,
                                  exclude_artifacts = FALSE,
                                  var_type = 1,
                                  backfit = TRUE, nullcheck = TRUE,
                                  covariates = NULL,
                                  covariate_mode = c("residualize", "fixed"),
                                  seed = 1, flash_control = list(),
                                  feature_scale = "unknown",
                                  preprocessing_spec = NULL,
                                  backend = c("flashier", "cebmf_torch"),
                                  var_floor = "none") {
  stopifnot(inherits(x, "cell_program_data"))
  covariate_mode <- match.arg(covariate_mode)
  backend <- match.arg(backend)
  # var_floor: "none" (default, behavior-preserving) | "vignette" (the flashier
  # single-cell vignette's fixed residual SD, S = sd(log(Poisson(1/n)+1)), n =
  # donors in the fit matrix; drawn once per cell type with a fixed local seed
  # so the value is deterministic and independent of the fit seed) | a numeric
  # scalar S. When supplied, the residual SDs are FIXED at S (flashier and
  # cebmf_torch both stop estimating variances) — see the CAND-13 demonstrator
  # for the calibration evidence (reports/audit/E_programs/).
  S_value <- NULL
  vignette_floor <- FALSE
  if (identical(var_floor, "vignette")) {
    vignette_floor <- TRUE
  } else if (is.numeric(var_floor)) {
    if (length(var_floor) != 1L || !is.finite(var_floor) || var_floor <= 0) {
      cli::cli_abort("{.arg var_floor} must be \"none\", \"vignette\", or a positive finite scalar.")
    }
    S_value <- as.numeric(var_floor)
  } else if (!identical(var_floor, "none")) {
    cli::cli_abort("{.arg var_floor} must be \"none\", \"vignette\", or a positive finite scalar.")
  }
  if (!is.null(exclude_genes) &&
      (!is.character(exclude_genes) || anyNA(exclude_genes))) {
    cli::cli_abort("{.arg exclude_genes} must be NULL or a character vector of gene IDs.")
  }
  if (isTRUE(exclude_artifacts)) {
    universe <- unique(unlist(lapply(x$pseudobulk, colnames), use.names = FALSE))
    exclude_genes <- union(exclude_genes, artifact_genes(universe))
  }
  if (covariate_mode == "fixed" && !is.null(covariates) && !center) {
    cli::cli_abort("{.code covariate_mode = \"fixed\"} requires {.code center = TRUE}.")
  }
  if (backend == "cebmf_torch" && covariate_mode == "fixed" && !is.null(covariates)) {
    cli::cli_abort("{.code covariate_mode = \"fixed\"} requires {.code backend = \"flashier\"}: the cebmf_torch backend has no fixed-loadings path. Use {.code covariate_mode = \"residualize\"} (the default) with this backend.")
  }
  if (!requireNamespace("flashier", quietly = TRUE) ||
      !requireNamespace("ebnm", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg flashier} and {.pkg ebnm} are required.")
  }
  # Ordinary fits route through flashier.utils::fit_ebmf() when available
  # (SYNC-1B); the fixed-covariate path (.flash_fixed) and degenerate sentinels
  # always use the legacy direct backend. Set
  # options(cellprograms.use_flashier_utils = FALSE) to force the legacy path.
  use_fu <- getOption("cellprograms.use_flashier_utils", TRUE) &&
    requireNamespace("flashier.utils", quietly = TRUE)
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

  # MR-0: a frozen preprocessing specification (freeze_preprocessing) governs
  # feature selection and adjustment. Conflicting explicit arguments are
  # refused (fail closed) so a refit can never silently drift from the frozen
  # panel; the spec digest is recorded in the fit's preprocessing provenance.
  frozen_digest <- NULL
  if (!is.null(preprocessing_spec)) {
    if (!inherits(preprocessing_spec, "cell_program_preproc_spec")) {
      cli::cli_abort("{.arg preprocessing_spec} must come from freeze_preprocessing().")
    }
    provided <- list(features = features, n_variable_genes = n_variable_genes,
      exclude_genes = exclude_genes, exclude_artifacts = exclude_artifacts,
      covariates = covariates, covariate_mode = covariate_mode,
      center = center, scale = scale, feature_scale = feature_scale)
    defaults <- list(features = "all", n_variable_genes = 2000, exclude_genes = NULL,
      exclude_artifacts = FALSE, covariates = NULL, covariate_mode = "residualize",
      center = TRUE, scale = FALSE, feature_scale = "unknown")
    conflicts <- names(provided)[vapply(names(provided),
      function(n) !identical(provided[[n]], defaults[[n]]), logical(1))]
    if (length(conflicts)) {
      cli::cli_abort("{.arg preprocessing_spec} is frozen: leave {conflicts} at their defaults; the spec governs.")
    }
    absent_ct <- setdiff(x$cell_types, names(preprocessing_spec$feature_sets))
    if (length(absent_ct)) {
      cli::cli_abort("Frozen spec lacks cell type(s): {.val {absent_ct}}.")
    }
    missing_genes <- setdiff(
      unique(unlist(preprocessing_spec$feature_sets)),
      unique(unlist(lapply(x$pseudobulk, colnames), use.names = FALSE)))
    if (length(missing_genes)) {
      cli::cli_abort("{length(missing_genes)} frozen gene(s) absent from this dataset, e.g. {.val {missing_genes[1]}}.")
    }
    feature_sets <- preprocessing_spec$feature_sets[x$cell_types]
    exclude_genes <- preprocessing_spec$exclude_genes
    covariates <- .spec_covariates(preprocessing_spec$covariates)
    covariate_mode <- preprocessing_spec$covariate_mode
    center <- preprocessing_spec$center
    scale <- preprocessing_spec$scale
    frozen_digest <- preprocessing_spec$digest
  }

  # R2/F3 (covariate audit): the fixed-mode centering guard above ran before
  # spec substitution; a spec frozen with covariate_mode = "fixed",
  # center = FALSE must not bypass it.
  if (covariate_mode == "fixed" && !is.null(covariates) && !center) {
    cli::cli_abort("{.code covariate_mode = \"fixed\"} requires {.code center = TRUE}.")
  }

  # Without a frozen spec, resolve the feature panel from the arguments.
  # With one, feature_sets was already taken from the spec above: re-running
  # .select_features() here would clobber a frozen "variable"-mode panel,
  # because the conflict check pins features/n_variable_genes at defaults.
  if (is.null(preprocessing_spec)) {
    feature_sets <- .select_features(x, features, n_variable_genes, exclude_genes)
  }

  fits <- list()
  for (ct in x$cell_types) {
    Y <- as.matrix(x$pseudobulk[[ct]])
    genes <- feature_sets[[ct]]
    Y <- Y[, genes, drop = FALSE]
    Xres <- if (!is.null(covariates) && covariate_mode == "residualize") {
      .cov_design(rownames(Y), x$sample_metadata, covariates, ct)
    } else NULL
    # MR-1 replay record: the exact adjustment applied to this cell type's
    # pseudobulk, so project_programs_to_cells() can put single-cell query
    # data into the same space (residualization coefficients, then the
    # centering/scale values actually used here).
    replay <- list(center = rep(0, ncol(Y)), scale = rep(1, ncol(Y)),
                   residualize_coefficients = NULL, design_columns = NULL)
    if (!is.null(Xres)) {
      res <- stats::lm.fit(Xres, Y)
      R <- res$residuals
      dimnames(R) <- dimnames(Y)
      Y <- R
      replay$residualize_coefficients <- res$coefficients
      replay$design_columns <- colnames(Xres)
    }
    if (center) {
      Y <- scale(Y, center = TRUE, scale = scale)
      replay$center <- unname(attr(Y, "scaled:center"))
      if (!is.null(attr(Y, "scaled:scale"))) replay$scale <- unname(attr(Y, "scaled:scale"))
    } else if (scale) {
      Y <- scale(Y, center = FALSE, scale = TRUE)
      replay$scale <- unname(attr(Y, "scaled:scale"))
    }
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
        flash = list(L_pm = NULL, F_pm = NULL, L_psd = NULL, F_psd = NULL),
        Y_used = Y,
        genes = character(0),
        elapsed_secs = 0,
        skipped = TRUE,
        preproc_replay = replay
      )
      next
    }
    if (!all(keep_cols)) {
      cli::cli_inform("Dropping {sum(!keep_cols)} zero-variance gene(s) in cell type {.val {ct}}.")
      Y <- Y[, keep_cols, drop = FALSE]
      genes <- genes[keep_cols]
      # R2/F1 (covariate audit): subset the replay vectors to the retained
      # genes so projection replay indexes the same positions as Y_used.
      replay$center <- replay$center[keep_cols]
      replay$scale  <- replay$scale[keep_cols]
    }

    # C4 safety guard (ROADMAP M6, D7c): the requested number of factors must
    # fit within the residual degrees of freedom. q_ct counts the covariate
    # columns adjusted for in this cell type (residualize: design rank minus
    # the intercept; fixed: the fixed columns), so n_ct - q_ct - 1 is the
    # degrees of freedom left after the covariate model and one spare df.
    # Never binds at operating points (n >= 50, max_factors <= 10); it fires
    # for small-n cell types at the CMP-detection K (e.g. max_factors = 30).
    q_ct <- if (!is.null(Xres)) {
      max(qr(Xres)$rank - 1L, 0L)
    } else if (!is.null(Xc)) {
      ncol(Xc)
    } else 0L
    n_ct <- nrow(Y)
    if (max_factors > n_ct - q_ct - 1L) {
      cli::cli_abort(c(
        "Too many factors in cell type {.val {ct}}: {.arg max_factors} = {max_factors} exceeds the {n_ct - q_ct - 1L} factor(s) the residual degrees of freedom support (n = {n_ct} donor(s), q = {q_ct} covariate column(s)).",
        i = "Lower {.arg max_factors}, drop covariates, or add donors for this cell type."
      ))
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

    # Per-cell-type vignette floor: n = donors in this cell type's fit matrix.
    # S_vf is the var_floor-derived S (distinct from a flash_control$S override,
    # which keeps the registered fail-closed rejection in the flashier.utils
    # fast path — only var_floor-derived S reroutes the fit).
    S_vf <- NULL
    if (vignette_floor) {
      set.seed(20261016)  # fixed local seed: S is deterministic, fit-seed-independent
      S_vf <- sd(log(stats::rpois(1e7, 1 / nrow(Y)) + 1))
    } else if (!is.null(S_value)) {
      S_vf <- S_value
    }
    if (!is.null(S_vf)) args$S <- S_vf

    t0 <- Sys.time()
    set.seed(seed)
    # Ordinary fits (no fixed covariates) go through flashier.utils::fit_ebmf()
    # when available; the argument mapping below is exact (verified to machine
    # precision against do.call(flashier::flash, args)). flash_control members
    # that fit_ebmf() does not expose are rejected explicitly rather than
    # silently dropped. The fixed-covariate path always uses .flash_fixed().
    # A var_floor-derived S is not mappable through fit_ebmf(): route the fit
    # through the legacy backend with a warning.
    ordinary <- is.null(Xc) && use_fu && is.null(S_vf)
    if (!is.null(S_vf) && use_fu && is.null(Xc)) {
      cli::cli_warn("{.arg var_floor} requires the legacy flashier backend for this fit (flashier.utils::fit_ebmf() does not expose S); routing through {.fn flashier::flash}.")
    }
    run_flash <- function(args) {
      if (!is.null(Xc)) return(.flash_fixed(args, Xc))
      if (ordinary) return(.fit_via_flashier_utils(args, seed, feature_scale))
      do.call(flashier::flash, args)
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
    if (backend == "cebmf_torch") {
      # M6 adoption: the cebmf_ext wrapper carries its own numerical guard set
      # (sebetahat floor, tau cap, seed-family fallback) and raises loudly on
      # objective divergence — no flashier-style fallback ladder here.
      fit <- tryCatch(.fit_via_cebmf(Y, max_factors, seed, S = args$S),
                      error = function(e) {
                        if (is.null(first_err)) first_err <<- e
                        NULL
                      })
      if (is.null(fit)) {
        cli::cli_abort("cebmf_torch failed in cell type {.val {ct}}: {conditionMessage(first_err)}")
      }
      if (!.cebmf_fit_ok(fit)) {
        cli::cli_abort("cebmf_torch returned a non-finite objective in cell type {.val {ct}}.")
      }
    } else {
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
    }
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!is.null(fit$fixed_L)) Y <- Y - fit$fixed_L %*% t(fit$fixed_F)
    fits[[ct]] <- list(
      flash = fit,
      Y_used = Y,
      genes = genes,
      elapsed_secs = elapsed,
      fallback = fallback,
      S_used = args$S,
      preproc_replay = replay
    )
  }

  x$fits <- fits
  x$preprocessing <- list(
    centered = center, scaled = scale,
    feature_sets = feature_sets,
    exclude_genes = exclude_genes,
    exclude_artifacts = isTRUE(exclude_artifacts),
    loading_prior = loading_prior,
    max_factors = max_factors, var_type = var_type,
    backfit = backfit, nullcheck = nullcheck,
    covariates = covariates, covariate_mode = covariate_mode,
    flash_control = flash_control,
    feature_scale = feature_scale,
    frozen_spec = frozen_digest,
    # CAND-12: store the FULL frozen spec, not just its digest — the digest
    # alone made replay impossible in principle and let stability refits
    # bypass the MR-0 fail-closed checks. NULL for unfrozen fits.
    spec = preprocessing_spec,
    backend = backend,
    var_floor = var_floor
  )
  x$provenance <- list(
    seed = seed,
    timestamp = Sys.time(),
    r_version = R.version.string,
    flashier_version = as.character(utils::packageVersion("flashier")),
    ebnm_version = as.character(utils::packageVersion("ebnm")),
    backend = backend,
    var_floor = var_floor
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
    # R2/F4 (covariate audit): drop every TERM that references a
    # within-ct-constant variable, including interaction terms — drop.terms()
    # on main-effect labels alone leaves e.g. age:sex in the formula and
    # model.matrix() crashes obscurely ("variable lengths differ").
    tl <- attr(stats::terms(f), "term.labels")
    tl <- tl[vapply(tl, function(t) !any(all.vars(stats::as.formula(paste("~", t))) %in% vars[const]),
                    logical(1))]
    f <- stats::reformulate(tl)   # all-constant case already returned NULL above
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
    L_psd = if (!is.null(fl$L_psd)) fl$L_psd[, k, drop = FALSE],
    F_psd = if (!is.null(fl$F_psd)) fl$F_psd[, k, drop = FALSE],
    F_lfsr = if (!is.null(fl$F_lfsr)) fl$F_lfsr[, k, drop = FALSE],
    fixed_L = Xc, fixed_F = fl$F_pm[, seq_len(q), drop = FALSE]
  )
  if (length(k) == 0L) { out$L_pm <- NULL; out$F_pm <- NULL }
  out
}

# Route an ordinary fit through flashier.utils::fit_ebmf(). `args` is the
# cellprograms flash argument list (data, ebnm_fn, var_type, greedy_Kmax,
# backfit, nullcheck, verbose, plus any flash_control overrides). The mapping
# to fit_ebmf() is exact; flash_control members that fit_ebmf() does not
# expose are rejected explicitly so no control is ever silently dropped.
# `feature_scale` is declared measurement units for the manifest (MR-1: the
# projection oracle refuses queries whose declared units differ from the
# source manifest, so fits must not stay "unknown" if they will be projected).
# It is passed separately rather than through `args` because the legacy
# backend do.call(flashier::flash, args) has no such formal.
.fit_via_flashier_utils <- function(args, seed, feature_scale = "unknown") {
  # flash_control overrides that fit_ebmf() cannot map: anything beyond the
  # core argument set. The core set is forwarded explicitly below.
  core <- c("data", "ebnm_fn", "var_type", "greedy_Kmax", "backfit",
            "nullcheck", "verbose")
  extra <- setdiff(names(args), core)
  if (length(extra)) {
    cli::cli_abort(c(
      "flash_control setting(s) {.val {extra}} cannot be mapped through flashier.utils::fit_ebmf().",
      i = "Set {.code options(cellprograms.use_flashier_utils = FALSE)} to use the legacy backend for this fit."
    ))
  }
  flashier.utils::fit_ebmf(
    X = args$data,
    sample_side = "rows",
    ebnm_fn = args$ebnm_fn,
    var_type = args$var_type,
    max_factors = args$greedy_Kmax,
    backfit = args$backfit,
    nullcheck = args$nullcheck,
    seed = seed,
    verbose = args$verbose,
    feature_scale = feature_scale
  )
}

# Feature selection: independent per cell type; never uses outcomes.
# `exclude_genes` is applied to the universe BEFORE selection, so excluded
# genes cannot consume variable-gene slots; with an explicit `features` list
# they are dropped from the requested vectors after validation.
.select_features <- function(x, features, n_variable_genes, exclude_genes = NULL) {
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
      setdiff(g, exclude_genes)
    }) |> stats::setNames(x$cell_types))
  }
  out <- list()
  for (ct in x$cell_types) {
    genes <- setdiff(colnames(x$pseudobulk[[ct]]), exclude_genes)
    if (identical(features, "all")) {
      out[[ct]] <- genes
    } else if (identical(features, "variable")) {
      v <- apply(as.matrix(x$pseudobulk[[ct]])[, genes, drop = FALSE], 2, stats::var)
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
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' head(fit$scores$B[, 1])     # canonicalized donor scores of program B_1
#' head(fit$loadings$B[, 1])   # canonicalized gene loadings of B_1
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
    # SYNC-2: ordinary fits carry a native flash object with flashbridge
    # metadata and route through flashier.utils::standardize_factors("rms"),
    # which reproduces the inline RMS + sign convention exactly (verified to
    # 0 difference). .flash_fixed() outputs and legacy objects lack the
    # metadata and use the inline branch below.
    use_upstream <- requireNamespace("flashier.utils", quietly = TRUE) &&
      !is.null(attr(f, "flashbridge_metadata"))
    if (use_upstream) {
      rms <- flashier.utils::standardize_factors(
        f, sample_side = "rows", scaling = "rms", orientation = "anchor_feature")
      Zc <- as.matrix(rms$activity)
      Wc <- as.matrix(rms$effects)
      # preserve the native dimnames so downstream reconstruction/extraction
      # matches the legacy inline branch exactly
      dimnames(Zc) <- dimnames(Z)
      dimnames(Wc) <- dimnames(W)
      s_k <- unname(rms$manifest$rms_scale)
      sign_k <- unname(rms$manifest$orientation_sign)
    } else {
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
    }
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
