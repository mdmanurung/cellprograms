# Maintained native-fit and actual-refit regressions, promoted from the bounded
# CORE-01/02/03 acceptance fixtures. Run by both test_local and installed checks.
# Tiny synthetic real-valued inputs are fixtures, not biological evidence.
local({
audit <- new.env(parent = emptyenv())
audit$evidence <- list()
audit$conditions <- list()

# Capture warnings/messages explicitly rather than hiding their existence.
observed <- function(label, expr) {
  ww <- mm <- character()
  err <- NULL
  value <- tryCatch(withCallingHandlers(
    force(expr),
    warning = function(w) { ww <<- c(ww, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { mm <<- c(mm, conditionMessage(m)); invokeRestart("muffleMessage") }),
    error = function(e) { err <<- e; NULL })
  audit$conditions[[label]] <- list(warnings = ww, messages = mm,
    error = if (is.null(err)) NULL else conditionMessage(err))
  list(value = value, error = err, warnings = ww, messages = mm)
}

toy <- function(n = 24L, p = 12L, seed = 3003L) {
  withr::local_seed(seed)
  age <- seq(20, 70, length.out = n)
  z <- sin(seq(0, 4 * pi, length.out = n))
  y <- outer(z, rep(c(5, -4, 3, -3), length.out = p)) +
    outer(as.numeric(scale(age)), rep(c(.4, -.2), length.out = p)) +
    matrix(rnorm(n * p, sd = .12), n, p)
  dimnames(y) <- list(paste0("d", seq_len(n)), paste0("g", seq_len(p)))
  as_cell_program_data(list(A = y), sample_metadata = data.frame(
    observation_id = rownames(y), age = age, rare = c(1, rep(0, n - 1L))))
}

tiny_fit <- function(x, ...) {
  do.call(fit_celltype_programs, c(list(x = x, max_factors = 1L,
    loading_prior = "point_normal", backfit = FALSE, nullcheck = FALSE,
    seed = 17L), list(...)))
}

expect_factors_equal <- function(a, b, tol = 1e-8) {
  expect_equal(dim(a$L_pm), dim(b$L_pm))
  expect_equal(dim(a$F_pm), dim(b$F_pm))
  expect_equal(as.numeric(a$L_pm), as.numeric(b$L_pm), tolerance = tol)
  expect_equal(as.numeric(a$F_pm), as.numeric(b$F_pm), tolerance = tol)
  # Kronecker noise is a list of row/column components, not a scalar/vector.
  expect_equal(a$residuals_sd, b$residuals_sd, tolerance = tol)
  expect_equal(a$elbo, b$elbo, tolerance = tol)
}

test_that("CORE-01 real ordinary configuration matrix agrees with native flashier", {
  x <- toy()
  rows <- list(
    list(id = "upstream-defaults", prior = "point_laplace", vt = 1, bf = TRUE, nc = TRUE),
    list(id = "shared-normal", prior = "point_normal", vt = 0, bf = FALSE, nc = FALSE),
    list(id = "column-unimodal", prior = "unimodal", vt = 2, bf = TRUE, nc = FALSE),
    list(id = "kronecker", prior = "point_normal", vt = c(1, 2), bf = FALSE, nc = TRUE),
    list(id = "known-fixed-binding", prior = "point_laplace", vt = 1, bf = FALSE, nc = FALSE,
         S = .4),
    list(id = "direct-estimated", prior = "point_normal", vt = 2, bf = TRUE, nc = TRUE,
         direct = TRUE),
    list(id = "fixed-covariate-known", prior = "point_normal", vt = 1, bf = FALSE,
         nc = FALSE, S = .35, mode = "fixed"),
    list(id = "residualized", prior = "point_normal", vt = 1, bf = FALSE,
         nc = FALSE, mode = "residualize")
  )
  for (r in rows) {
    options(cellprograms.use_flashier_utils = !isTRUE(r$direct))
    on.exit(options(cellprograms.use_flashier_utils = NULL), add = TRUE)
    args <- list(x = x, max_factors = 1L, loading_prior = r$prior,
      var_type = r$vt, backfit = r$bf, nullcheck = r$nc, seed = 17L,
      feature_scale = "synthetic_linear", var_floor = if (is.null(r$S)) "none" else r$S,
      flash_control = list(verbose = 0L))
    if (!is.null(r$mode)) {
      args$covariates <- ~age
      args$covariate_mode <- r$mode
    }
    got <- observed(r$id, do.call(fit_celltype_programs, args))
    expect_null(got$error, info = r$id)
    if (!is.null(got$error)) next
    f <- got$value$fits$A
    e <- f$effective_config
    expect_identical(f$fallback, "none", info = r$id)
    expect_identical(e$backend, "flashier")
    expect_identical(e$loading_prior, r$prior)
    expect_identical(e$sample_prior, "normal")
    expect_identical(e$backfit, r$bf)
    expect_identical(e$nullcheck, r$nc)
    expect_equal(e$max_factors, 1)
    expect_equal(e$var_type, if (is.null(r$S)) r$vt else NULL)
    expect_equal(e$noise_sd, r$S)
    expect_identical(e$noise_mode, if (is.null(r$S)) "estimated" else "known_fixed")
    engine <- if (identical(r$mode, "fixed")) "flashier_fixed" else if (
      isTRUE(r$direct) || !is.null(r$S)) "flashier" else "flashier.utils"
    expect_identical(e$engine_path, engine)
    raw <- x$pseudobulk$A
    if (identical(r$mode, "residualize")) {
      raw <- stats::lm.fit(model.matrix(~age, x$sample_metadata), raw)$residuals
      dimnames(raw) <- dimnames(x$pseudobulk$A)
    }
    y <- scale(raw, center = TRUE, scale = FALSE)
    prior <- switch(r$prior, point_normal = ebnm::ebnm_point_normal,
                     point_laplace = ebnm::ebnm_point_laplace,
                     unimodal = ebnm::ebnm_unimodal)
    native <- observed(paste0(r$id, "-native"), {
      withr::local_seed(17L)
      if (identical(r$mode, "fixed")) {
        xc <- matrix(as.numeric(scale(x$sample_metadata$age)), ncol = 1L)
        f0 <- t(solve(crossprod(xc), crossprod(xc, y)))
        fl <- flashier::flash_init(y, S = r$S, var_type = NULL)
        fl <- flashier::flash_set_verbose(fl, 0L)
        fl <- flashier::flash_factors_init(fl, list(xc, f0),
          ebnm_fn = flashier::flash_ebnm(prior_family = "normal", mode = 0, scale = 10))
        fl <- flashier::flash_factors_fix(fl, kset = 1L, which_dim = "loadings")
        fl <- flashier::flash_greedy(fl, Kmax = 1L,
                                     ebnm_fn = list(ebnm::ebnm_normal, prior))
        free <- setdiff(seq_len(fl$n_factors), 1L)
        list(L_pm = fl$L_pm[, free, drop = FALSE], F_pm = fl$F_pm[, free, drop = FALSE],
          residuals_sd = fl$residuals_sd, elbo = fl$elbo)
      } else {
        flashier::flash(y, S = r$S, var_type = if (is.null(r$S)) r$vt else NULL,
          ebnm_fn = list(ebnm::ebnm_normal, prior), greedy_Kmax = 1L,
          backfit = r$bf, nullcheck = r$nc, verbose = 0L)
      }
    })
    expect_null(native$error, info = r$id)
    if (!is.null(native$error)) next
    expect_factors_equal(f$flash, native$value)
    if (!is.null(r$S)) {
      expect_equal(as.numeric(f$flash$residuals_sd), r$S, tolerance = 1e-12)
      expect_equal(f$S_used, r$S)
    }
    if (engine == "flashier.utils") {
      manifest <- attr(f$flash, "flashbridge_metadata")
      expect_equal(manifest$configuration$var_type, e$var_type)
      expect_identical(manifest$configuration$backfit, e$backfit)
      expect_identical(manifest$configuration$nullcheck, e$nullcheck)
      expect_identical(manifest$feature_scale, "synthetic_linear")
      expect_identical(manifest$prior_identifiers[2L], paste0("ebnm_", e$loading_prior))
    }
    audit$evidence[[r$id]] <- list(requested = r, effective = e,
      native_elbo = native$value$elbo, actual_elbo = f$flash$elbo,
      residuals_sd = f$flash$residuals_sd, fallback = f$fallback)
  }
})

test_that("CORE-01 unsafe controls and optional estimator defaults refuse before any fit", {
  x <- toy()
  # This spy must never be reached; unlike successful fixtures it performs no fitting.
  local_mocked_bindings(
    .fit_via_flashier_utils = function(...) stop("ENGINE REACHED"),
    .flash_fixed = function(...) stop("ENGINE REACHED"),
    .fit_via_cebmf = function(...) stop("ENGINE REACHED"), .package = "cellprograms")
  for (control in list(list(backfit = FALSE), list(nullcheck = FALSE),
    list(var_type = 2), list(S = .2), list(greedy_Kmax = 2),
    list(ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_normal)), list(0L))) {
    expect_error(tiny_fit(x, flash_control = control), "cannot be mapped safely")
  }
  expect_error(fit_celltype_programs(x, backend = "cebmf_torch"),
               "distinct experimental")
  expect_error(fit_celltype_programs(x, loading_prior = "laplace"),
               "only supported")
  expect_error(fit_celltype_programs(x, max_factors = 1, backend = "cebmf_torch",
    loading_prior = "laplace", var_type = 2, nullcheck = FALSE,
    covariates = ~age, covariate_mode = "fixed"), "no fixed-loadings")
})

test_that("CORE-02 old schema roundtrip refits unknown units and preserves gene identity", {
  x <- toy()
  s <- freeze_preprocessing(x, features = list(A = rev(colnames(x$pseudobulk$A))),
                            feature_scale = "synthetic_linear")
  old <- s
  old$digest <- NULL
  old$spec_version <- "1.0.0"
  old$feature_scale <- NULL
  old$digest <- digest::digest(old, algo = "sha256")
  old <- unserialize(serialize(old, NULL, version = 3))
  f <- tiny_fit(x, preprocessing_spec = old)
  expect_identical(f$preprocessing$feature_scale, "unknown")
  expect_identical(f$preprocessing$spec$spec_version, "1.0.0")
  expect_identical(f$fits$A$genes, old$feature_sets$A)
  # Query matrix storage order must not change the frozen fit's feature identities.
  xr <- x
  xr$pseudobulk$A <- xr$pseudobulk$A[, rev(colnames(xr$pseudobulk$A)), drop = FALSE]
  fr <- tiny_fit(xr, preprocessing_spec = old)
  expect_factors_equal(fr$fits$A$flash, f$fits$A$flash)
  stale <- old
  stale$feature_sets$A <- rev(stale$feature_sets$A)
  expect_error(tiny_fit(x, preprocessing_spec = stale), "digest mismatch")
  future <- old
  future$digest <- NULL
  future$spec_version <- "2.0.0"
  future$digest <- digest::digest(future, algo = "sha256")
  expect_error(tiny_fit(x, preprocessing_spec = future), "Unsupported")
  audit$evidence$old_schema <- list(input_schema = old$spec_version,
    output_units = f$preprocessing$feature_scale, output_genes = f$fits$A$genes,
    note = "Synthetic v1.0 shape, not an archived historical binary or migration certificate.")
})

test_that("CORE-02 complete and legacy linear/transformed fits serialize with honest refusal", {
  x <- toy()
  for (form in list(~age, ~poly(age, 2))) {
    label <- paste(deparse(form), collapse = "")
    fit <- canonicalize_programs(tiny_fit(x, covariates = form,
                                          feature_scale = "synthetic_linear"))
    copy <- unserialize(serialize(fit, NULL, version = 3))
    expect_equal(copy$scores, fit$scores, tolerance = 0)
    expect_equal(copy$loadings, fit$loadings, tolerance = 0)
    idx <- c(7L, 2L, 7L)
    q <- x$pseudobulk$A[idx, , drop = FALSE]
    rownames(q) <- paste0("q", seq_along(idx))
    ct <- setNames(rep("A", length(idx)), rownames(q))
    donors <- setNames(rownames(x$pseudobulk$A)[idx], rownames(q))
    args <- list(expression = q, cell_types = ct, donors = donors,
      sample_metadata = x$sample_metadata)
    pr <- observed(paste0("current-projection-", label),
                    do.call(project_programs_to_cells, c(list(fit = copy), args)))
    expect_null(pr$error)
    if (!is.null(pr$error)) next
    expect_true(all(is.finite(pr$value$activity)))
    expect_equal(pr$value$activity[1, ], pr$value$activity[3, ], tolerance = 1e-10)
    # Deleting the trained state emulates an old fit shape. Both linear and
    # transformed legacy residualized objects require refit under current code.
    legacy <- copy
    legacy$fits$A$preproc_replay$design_state <- NULL
    legacy <- unserialize(serialize(legacy, NULL, version = 3))
    expect_equal(legacy$scores, fit$scores, tolerance = 0)
    expect_error(do.call(project_programs_to_cells, c(list(fit = legacy), args)),
                 "refit legacy residualized models")
    missing <- copy
    missing$fits$A$preproc_replay <- NULL
    expect_error(do.call(project_programs_to_cells, c(list(fit = missing), args)),
                 "predates preprocessing replay records")
    legacy_backend <- copy
    attr(legacy_backend$fits$A$flash, "flashbridge_metadata") <- NULL
    expect_error(do.call(project_programs_to_cells,
                        c(list(fit = legacy_backend), args)), "no flashier.utils manifest")
    expect_error(do.call(project_programs_to_cells,
      c(list(fit = copy, feature_scale = "different_synthetic_units"), args)),
      "differs from the fit's[[:space:]]+recorded units")
    audit$evidence[[paste0("legacy-", label)]] <- list(
      complete_roundtrip_projection = "finite; repeated donor invariant",
      missing_state_projection = "refit required, including linear formula",
      missing_replay_projection = "refit required", no_manifest_projection = "refit required",
      old_factors_unchanged = TRUE)
  }
})

test_that("CORE-02 TR03-D001 unknown references require a declared-unit refit", {
  x <- toy()
  f <- canonicalize_programs(tiny_fit(x)) # no units inferred from values
  q <- x$pseudobulk$A[1:3, , drop = FALSE]
  rownames(q) <- paste0("q", 1:3)
  ct <- setNames(rep("A", 3), rownames(q))
  unchanged <- serialize(f, NULL)
  expect_error(project_programs_to_cells(f, q, ct), "measurement units are undeclared")
  # An explicit owner-delegated scope amendment, not a retroactive pass of the
  # old advertised override. Neither caller nor implementation relabels f.
  attempt <- observed("explicit-unknown-units",
    project_programs_to_cells(f, q, ct, feature_scale = "synthetic_linear"))
  audit$evidence$unknown_units <- list(
    fit_units = f$preprocessing$feature_scale,
    native_units = attr(f$fits$A$flash, "flashbridge_metadata")$feature_scale,
    explicit_truthful_query_units = "synthetic_linear",
    error = if (is.null(attempt$error)) NULL else conditionMessage(attempt$error))
  expect_s3_class(attempt$error, "error")
  expect_match(conditionMessage(attempt$error), "refit")
  expect_identical(serialize(f, NULL), unchanged)
  declared <- canonicalize_programs(tiny_fit(x, feature_scale = "synthetic_linear"))
  projection <- project_programs_to_cells(declared, q, ct,
                                           feature_scale = "synthetic_linear")
  expect_true(all(is.finite(projection$activity)))
})

test_that("CORE-03 frozen and unfrozen stability actually replay the same ordinary estimator", {
  x <- toy()
  fit_original <- get("fit_celltype_programs", envir = asNamespace("cellprograms"))
  for (frozen in c(FALSE, TRUE)) {
    args <- list(feature_scale = "synthetic_linear", covariates = ~age)
    if (frozen) {
      spec <- freeze_preprocessing(x, features = list(A = rev(colnames(x$pseudobulk$A))),
        feature_scale = "synthetic_linear", covariates = ~age)
      args <- list(preprocessing_spec = spec)
    }
    reference <- observed(paste0("stability-reference-", frozen),
      canonicalize_programs(do.call(tiny_fit, c(list(x = x, var_floor = .4), args))))
    expect_null(reference$error)
    if (!is.null(reference$error)) next
    ref <- reference$value
    captures <- list()
    wrapper <- function(...) {
      aa <- list(...)
      ans <- do.call(fit_original, aa)
      captures[[length(captures) + 1L]] <<- list(args = aa, fit = ans)
      ans
    }
    # A transparent instrumented REAL fit, not a canned fit or forwarding-only mock.
    run <- function() {
      local_mocked_bindings(fit_celltype_programs = wrapper, .package = "cellprograms")
      assess_program_stability(ref, n_boot = 2L, resample_frac = .75,
                                seed = 31L, progress = FALSE)
    }
    got <- observed(paste0("stability-identity-", frozen), run())
    expect_null(got$error)
    if (!is.null(got$error)) next
    st <- got$value$stability
    expect_identical(st$valid_refits, 2L)
    expect_identical(st$failed_refits, 0L)
    expect_length(captures, 2L)
    expect_true(all(st$replicate_log$status == "ok"))
    actual <- list()
    for (b in seq_along(captures)) {
      cf <- captures[[b]]$fit
      direct <- observed(paste0("independent-replay-", frozen, "-", b),
                          do.call(fit_original, captures[[b]]$args))
      expect_null(direct$error)
      if (!is.null(direct$error)) next
      expect_factors_equal(cf$fits$A$flash, direct$value$fits$A$flash, tol = 0)
      expect_identical(cf$preprocessing$feature_scale, "synthetic_linear")
      expect_identical(cf$fits$A$effective_config$engine_path, "flashier")
      expect_identical(cf$fits$A$effective_config$loading_prior, "point_normal")
      expect_identical(cf$fits$A$effective_config$noise_mode, "known_fixed")
      expect_identical(cf$fits$A$fallback, "none")
      expect_equal(as.numeric(cf$fits$A$flash$residuals_sd), .4, tolerance = 1e-12)
      expect_identical(st$replicate_log$effective_refit[[b]]$A$effective_config,
                       cf$fits$A$effective_config)
      expect_equal(st$replicate_log$effective_refit[[b]]$A$S_used, .4)
      expect_equal(st$replicate_log$n_obs[[b]][["A"]], 18)
      expect_equal(cf$provenance$seed, 17L + b)
      if (frozen) expect_identical(cf$preprocessing$frozen_spec, spec$digest)
      if (!frozen) expect_null(cf$preprocessing$frozen_spec)
      actual[[b]] <- list(row_ids = rownames(cf$pseudobulk$A),
        effective = cf$fits$A$effective_config, units = cf$preprocessing$feature_scale,
        frozen_digest = cf$preprocessing$frozen_spec, fallback = cf$fits$A$fallback,
        residuals_sd = cf$fits$A$flash$residuals_sd)
    }
    expect_equal(st$per_program$n_eligible, rep(2, nrow(st$per_program)))
    expect_equal(st$per_program$recovery_freq, st$per_program$n_hit / 2)
    audit$evidence[[paste0("stability-identity-", frozen)]] <- list(
      actual = actual, replicate_log = st$replicate_log, per_program = st$per_program,
      note = "Run externally source-bound by receipt; native fit provenance lacks a source digest.")
  }
})

test_that("CORE-03 genuine partial failures retain requested successful and eligible denominators", {
  x <- toy(n = 12L)
  ref <- observed("mixed-reference", canonicalize_programs(
    tiny_fit(x, covariates = ~rare, feature_scale = "synthetic_linear")))
  expect_null(ref$error)
  if (!is.null(ref$error)) return(invisible(NULL))
  run <- observed("mixed-real-refits", assess_program_stability(ref$value, n_boot = 6L,
    resample_frac = .75, seed = 31L, progress = FALSE))
  expect_null(run$error)
  if (!is.null(run$error)) return(invisible(NULL))
  st <- run$value$stability
  log <- st$replicate_log
  # Actual selected rows determine whether the rare covariate is constant.
  # Five failures are real pre-fit residual-df refusals, not injected errors.
  expect_identical(log$status, c("error", "error", "error", "error", "ok", "error"))
  expect_identical(st$valid_refits, 1L)
  expect_identical(st$failed_refits, 5L)
  expect_equal(st$valid_refits + st$failed_refits, st$n_boot)
  expect_equal(nrow(log), st$n_boot)
  expect_true(all(grepl("need >= 10 residual", log$error_message[log$status == "error"])))
  expect_true(any(grepl("5 of 6 subsample refits failed", run$warnings)))
  expect_true(all(vapply(log$effective_refit[log$status == "error"], is.null, logical(1))))
  pp <- st$per_program
  expect_equal(pp$n_eligible, rep(1, nrow(pp)))
  expect_equal(pp$recovery_freq, pp$n_hit / 6)
  expect_equal(pp$conditional_recovery_freq, pp$n_hit / st$valid_refits)
  audit$evidence$mixed_denominators <- list(log = log, per_program = pp,
    requested = st$n_boot, valid = st$valid_refits, failed = st$failed_refits)
})

test_that("CORE-03 successful K-zero refits are not failures and not program eligible", {
  x <- toy(n = 12L)
  x$pseudobulk$A[,] <- 0
  x$pseudobulk$A[1L, ] <- rep(c(10, -10, 5, -5), 3L)
  reference <- observed("zero-factor-reference",
    canonicalize_programs(tiny_fit(x, var_floor = .4, feature_scale = "synthetic_linear")))
  expect_null(reference$error)
  if (!is.null(reference$error)) return(invisible(NULL))
  ref <- reference$value
  expect_equal(ncol(ref$scores$A), 1)
  run <- observed("zero-factor-real-refits", assess_program_stability(ref, n_boot = 6L,
    resample_frac = .75, seed = 31L, progress = FALSE))
  expect_null(run$error)
  if (!is.null(run$error)) return(invisible(NULL))
  st <- run$value$stability
  k <- vapply(st$replicate_log$n_factors, `[[`, integer(1), "A")
  expect_identical(st$valid_refits, 6L)
  expect_identical(st$failed_refits, 0L)
  expect_identical(k, c(1L, 1L, 1L, 1L, 0L, 1L))
  expect_identical(st$replicate_log$eligible, k > 0L)
  expect_equal(st$per_program$n_eligible, 5)
  expect_equal(st$per_program$n_hit, 5)
  expect_equal(st$per_program$recovery_freq, 5 / 6)
  expect_equal(st$per_program$conditional_recovery_freq, 5 / 6)
  audit$evidence$zero_factor_denominators <- list(log = st$replicate_log,
    per_program = st$per_program, requested = 6L, valid = 6L, failed = 0L)
})

test_that("CORE-03 actual corrupted-frozen-spec failures yield NA not zero stability", {
  x <- toy()
  spec <- freeze_preprocessing(x, feature_scale = "synthetic_linear")
  ref <- canonicalize_programs(tiny_fit(x, preprocessing_spec = spec))
  ref$preprocessing$spec$feature_sets$A <- rev(ref$preprocessing$spec$feature_sets$A)
  run <- observed("corrupt-frozen-real-refits", assess_program_stability(ref, n_boot = 2L,
    resample_frac = .75, seed = 31L, progress = FALSE))
  expect_null(run$error)
  if (!is.null(run$error)) return(invisible(NULL))
  st <- run$value$stability
  expect_identical(st$valid_refits, 0L)
  expect_identical(st$failed_refits, 2L)
  expect_equal(nrow(st$replicate_log), 2)
  expect_true(all(st$replicate_log$status == "error"))
  expect_true(all(grepl("digest mismatch", st$replicate_log$error_message)))
  expect_true(all(is.na(st$per_program$recovery_freq)))
  expect_true(all(is.na(st$per_program$conditional_recovery_freq)))
  expect_true(any(grepl("all 2 subsample refits failed", run$warnings)))
  expect_equal(st$per_program$n_eligible, 0)
  audit$evidence$all_fail <- list(log = st$replicate_log, per_program = st$per_program,
    requested = 2L, valid = 0L, failed = 2L)
})
})
