# cebmf_torch backend bridge (M6 adoption; ledger 2026-10-18).
#
# Layering rule (binding): R code here consumes ONLY the cebmf_ext result
# surface (the FitResult dataclass returned by fit_cebmf), never cebmf_torch
# internals — python/src/cebmf_ext is the only permitted entry point, and the
# copy shipped under inst/python/src is what production loads.
#
# Python resolution order:
#   1. option("cellprograms.cebmf_python_env"): a python binary or a venv dir
#      with cebmf_ext importable (used by the audit test suite).
#   2. basilisk (if installed): the pinned `cellprograms_cebmf` env
#      (torch==2.7.1 CPU + cebmf_torch@34abac5; the torch pin is load-bearing —
#      2.14.1+cpu diverges at iteration 0 on inputs where 2.7.1+cpu is green,
#      ledger 2026-10-18 P2 finding).
#   3. informative abort (fail closed; no silent fallback to flashier).
#
# Result mapping contract (registered, M6 review):
#   - L_pm / F_pm: subset to non-degenerate columns (column norm > 1e-8 * max).
#   - L_psd / F_psd / F_lfsr: NULL — the backend does not produce them and the
#     package never manufactures uncertainty it did not compute.
#   - elbo: NULL — cebmf_torch records a LOSS trace (negative ELBO + KL),
#     not an ELBO; surfacing it as one would mislabel the objective.
#     program_summary() shows NA via the %||% fallback; .fit_ok() passes on NULL.
#   - objective_trace / objective_direction / effective_k / converged / tau /
#     guard_diagnostics: carried through for provenance and gating.

.cebmf_backend_name <- "cebmf_torch"

.cebmf_basilisk_env <- function() {
  if (!requireNamespace("basilisk", quietly = TRUE)) {
    cli::cli_abort("{.pkg basilisk} is not installed.")
  }
  basilisk::BasiliskEnvironment(
    envname = "cellprograms_cebmf",
    pkgname = "cellprograms",
    packages = c("python=3.11"),
    pip = c(
      # The CPU-only torch wheel (+cpu local version) lives on the pytorch
      # wheel index, not PyPI; the index is supplied via PIP_EXTRA_INDEX_URL
      # in .cebmf_basilisk_start(). It cannot travel as a pip entry: basilisk
      # 1.18 passes pip entries through system2 unquoted (the shell splits
      # them) and basilisk 1.24 shell-quotes each entry (reticulate
      # virtualenv_create), so a flag entry breaks one or the other. An
      # environment variable passes through both paths. PyPI stays reachable
      # for the remaining entries.
      "torch==2.7.1+cpu",
      # cebmf_ext imports numpy and scipy directly (io.py, project.py);
      # torch does not pull either, so they must be pinned here. pandas is
      # the same story: cebmf_torch's fit path imports it but does not
      # declare it (surface: ModuleNotFoundError mid-fit, not at import).
      "numpy==2.4.6",
      "scipy==1.17.1",
      "pandas==3.0.6",
      # PEP 508 direct reference, commit-pinned (same pin as the venv path).
      # No spaces around "@": basilisk passes pip entries through system2
      # unquoted, so the shell would split a spaced reference into
      # "cebmf_torch", "@", "git+..." and pip would reject the bare "@".
      paste0(
        "cebmf_torch@git+https://github.com/william-denault/cebmf_torch.git@",
        "34abac5a7d05b6fc5c6f8052c5fc4fd96628597d"
      )
    )
  )
}

.cebmf_basilisk_start <- function(env) {
  # basilisk's pip-spec check greps every pip entry for "==" and cannot parse
  # PEP 508 direct references, so the commit-pinned cebmf_torch entry is
  # rejected even though a commit pin is fully reproducible. Relax the check
  # only while the environment is created (one-time; basilisk re-checks
  # nothing once the env directory exists). Upstream marks this knob
  # development-phase-only; the deviation is deliberate and documented in
  # reports/audit/TASKS.md.
  if (basilisk::getBasiliskCheckVersions()) {
    on.exit(basilisk::setBasiliskCheckVersions(TRUE), add = TRUE)
    basilisk::setBasiliskCheckVersions(FALSE)
  }
  # The CPU torch wheel needs the pytorch wheel index. A pip entry cannot
  # carry the index flag (basilisk 1.18: unquoted system2, shell splitting;
  # basilisk 1.24: reticulate shell-quotes each entry) - pip reads the index
  # from the environment instead, which passes through both code paths.
  if (!nzchar(Sys.getenv("PIP_EXTRA_INDEX_URL", ""))) {
    on.exit(Sys.unsetenv("PIP_EXTRA_INDEX_URL"), add = TRUE)
    Sys.setenv(PIP_EXTRA_INDEX_URL = "https://download.pytorch.org/whl/cpu")
  }
  basilisk::basiliskStart(env)
}

.cebmf_use_python_env <- function(opt) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_abort("The cebmf_torch backend requires the {.pkg reticulate} package.")
  }
  py <- if (file.exists(file.path(opt, "bin", "python"))) {
    file.path(opt, "bin", "python")
  } else {
    opt
  }
  if (!file.exists(py)) {
    cli::cli_abort("{.option cellprograms.cebmf_python_env} points to a non-existent python: {.val {opt}}.")
  }
  reticulate::use_python(py, required = TRUE)
  invisible(TRUE)
}

.cebmf_src_path <- function() {
  src <- system.file("python", "src", package = "cellprograms")
  if (!nzchar(src) || !dir.exists(file.path(src, "cebmf_ext"))) {
    cli::cli_abort("cebmf_ext source not found in the installed package (inst/python/src).")
  }
  src
}

.cebmf_fit_impl <- function(Y, max_factors, seed, S) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_abort("The cebmf_torch backend requires the {.pkg reticulate} package.")
  }
  reticulate::py_run_string(
    sprintf("import sys; sys.path.insert(0, r'%s')", .cebmf_src_path())
  )
  mod <- reticulate::import("cebmf_ext", convert = FALSE)
  torch <- reticulate::import("torch", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  Yt <- torch$from_numpy(np$array(Y)$astype("float32"))
  cfg <- mod$FitConfig(
    K = as.integer(max_factors),
    prior_L = "laplace", prior_F = "laplace",
    noise = "column_wise",
    allow_backfitting = TRUE,
    seed = as.integer(seed),
    S = S  # NULL -> None (noise estimated); scalar -> fixed KNOWN-noise SD
  )
  res <- mod$fit_cebmf(Yt, cfg)

  Lr <- reticulate::py_to_r(res$L)
  Fr <- reticulate::py_to_r(res$F)
  tau <- reticulate::py_to_r(res$tau)
  objective <- as.numeric(reticulate::py_to_r(res$objective))
  direction <- reticulate::py_to_r(res$direction)
  eff_k <- as.integer(reticulate::py_to_r(res$effective_k))
  converged <- reticulate::py_to_r(res$converged)
  gdiag <- reticulate::py_to_r(res$guard_diagnostics)

  # Non-degenerate column subset (registered mapping contract).
  norms <- sqrt(colSums(Lr * Lr))
  keep <- norms > 1e-8 * max(max(norms), 1e-300)
  K <- sum(keep)

  list(
    n_factors = K,
    elbo = NULL,                       # loss trace is not an ELBO (contract)
    objective_trace = objective,
    objective_direction = direction,
    residuals_sd = NULL,
    L_pm = Lr[, keep, drop = FALSE],
    F_pm = Fr[, keep, drop = FALSE],
    L_psd = NULL, F_psd = NULL, F_lfsr = NULL,  # never manufacture (contract)
    fixed_L = NULL, fixed_F = NULL,
    tau = tau,
    effective_k = eff_k,
    converged = converged,
    guard_diagnostics = gdiag,
    backend = .cebmf_backend_name
  )
}

.fit_via_cebmf <- function(Y, max_factors, seed, S = NULL) {
  stopifnot(is.matrix(Y))
  opt <- getOption("cellprograms.cebmf_python_env")
  if (!is.null(opt) && nzchar(opt)) {
    .cebmf_use_python_env(opt)
    .cebmf_fit_impl(Y, max_factors, seed, S)
  } else if (requireNamespace("basilisk", quietly = TRUE)) {
    env <- .cebmf_basilisk_env()
    cl <- .cebmf_basilisk_start(env)   # camelCase: basilisk exports basiliskStart/Run/Stop
    on.exit(basilisk::basiliskStop(cl), add = TRUE)
    basilisk::basiliskRun(cl, fun = .cebmf_fit_impl,
                          Y = Y, max_factors = max_factors, seed = seed, S = S)
  } else {
    cli::cli_abort(
      "No python environment for the cebmf_torch backend. Set {.option cellprograms.cebmf_python_env} to a python binary or venv with {.pkg cebmf_ext} importable, or install {.pkg basilisk} for the pinned environment."
    )
  }
}

.cebmf_fit_ok <- function(fit) {
  tr <- fit$objective_trace
  if (is.null(tr) || !length(tr)) return(FALSE)
  all(is.finite(tr)) && is.finite(sum(fit$L_pm)) && is.finite(sum(fit$F_pm))
}
