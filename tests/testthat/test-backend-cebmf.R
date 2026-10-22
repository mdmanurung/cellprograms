# Tests for the cebmf_torch backend bridge (M6 adoption), the var_floor
# exposure, and the R2 covariate-audit patches (F1-F4).
#
# Reconstructed 2026-10-18 (ledger): the original file was lost with the
# worker-1 repo tree during a sync; rewritten against the registered contract.
# Backend tests need a python env: option("cellprograms.cebmf_python_env")
# (set by audit runners) or basilisk; they skip otherwise. Covariate-patch and
# var_floor-validation tests are pure R.

cebmf_env_ready <- function() {
  opt <- getOption("cellprograms.cebmf_python_env")
  (!is.null(opt) && nzchar(opt)) || requireNamespace("basilisk", quietly = TRUE)
}

# Registered fixture recipe (replicates python/tests/conftest.py make_counts
# + ShiftedLog): block-structured counts with planted programs, shifted-log
# normalized and gene z-scaled. Raw-scale fixtures are OUT of the package's
# input regime and diverge the column-wise noise estimation (observed
# 2026-10-18; do not replace this recipe with unnormalized rnorm data).
.synth_counts <- function(n = 40, g = 50, p = 2, seed = 3) {
  set.seed(seed)
  genes_per_program <- max(3L, g %/% (2L * p))
  W <- matrix(0, g, p)
  program_genes <- list()
  for (j in seq_len(p)) {
    idx <- sample.int(g, genes_per_program)
    W[idx, j] <- stats::runif(genes_per_program, 0.5, 1.5)
    program_genes[[j]] <- idx
  }
  S <- matrix(0, n, p)
  for (j in seq_len(p)) {
    sel <- sample.int(n, n %/% 3)
    S[sel, j] <- stats::runif(length(sel), 1.0, 2.5)
  }
  lam <- exp(pmin(pmax(S %*% t(W), 0), 4)) + 0.05
  counts <- matrix(stats::rpois(length(lam), as.vector(lam)), n, g)
  sf <- rowSums(counts); sf <- sf / mean(sf)
  Y <- log1p(counts / sf)
  Y <- scale(Y)
  dimnames(Y) <- list(paste0("s", seq_len(n)), paste0("g", seq_len(g)))
  list(Y = Y, program_genes = program_genes, S_true = S)
}

# ---- var_floor validation (pure R) -----------------------------------------

test_that("var_floor rejects invalid values with the registered message", {
  g <- 20
  Y <- matrix(rnorm(10 * g), 10, g, dimnames = list(paste0("s", 1:10), paste0("g", 1:g)))
  cpd <- as_cell_program_data(list(CT1 = Y))
  for (bad in list("bogus", -1, c(0.5, 0.5), NA_real_)) {
    expect_error(fit_celltype_programs(cpd, var_floor = bad),
                 "must be \"none\", \"vignette\", or a positive finite scalar",
                 fixed = TRUE)
  }
})

# ---- R2/F3: fixed-mode centering guard after spec substitution --------------

test_that("F3: a frozen spec with covariate_mode=fixed, center=FALSE aborts", {
  set.seed(4)
  n <- 30; g <- 25
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  meta <- data.frame(observation_id = paste0("s", 1:n), age = seq(20, 60, length.out = n))
  cpd <- as_cell_program_data(list(CT1 = Y), sample_metadata = meta)
  # freeze_preprocessing itself does not guard the combination; the refit must.
  spec <- freeze_preprocessing(cpd, covariates = ~age, covariate_mode = "fixed", center = FALSE)
  expect_error(fit_celltype_programs(cpd, preprocessing_spec = spec),
               "center = TRUE", fixed = TRUE)
})

# ---- R2/F4: interaction terms referencing constant covariates are dropped ---

test_that("F4: interaction terms referencing a within-ct-constant covariate are dropped", {
  meta <- data.frame(
    observation_id = paste0("d", 1:20),
    age = seq(20, 60, length.out = 20),
    sex = factor(rep(c("F", "M"), 10)),
    site = factor(rep("A", 20))   # constant within this cell type
  )
  expect_message(
    X <- .cov_design(paste0("d", 1:20), meta, ~ age + sex + site + age:site, "ct"),
    "constant"
  )
  expect_false(any(grepl("site", colnames(X))))
  expect_true(all(c("(Intercept)", "age") %in% colnames(X)))
})

# ---- R2/F1: gene-drop block subsets the replay vectors ----------------------

test_that("F1: zero-variance gene drop subsets replay center/scale to retained genes", {
  set.seed(5)
  n <- 30; g <- 40
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  Y <- cbind(Y, gconst = 5)   # constant gene: centered away, then dropped
  cpd <- as_cell_program_data(list(CT1 = Y))
  fit <- fit_celltype_programs(cpd, max_factors = 3, seed = 1)
  ct_fit <- fit$fits$CT1
  expect_false("gconst" %in% ct_fit$genes)
  expect_equal(length(ct_fit$preproc_replay$center), length(ct_fit$genes))
  expect_equal(length(ct_fit$preproc_replay$scale), length(ct_fit$genes))
})

# ---- R2/F2: NA coefficients from a rank-deficient replay design are zeroed --

test_that("F2: projection with aliased covariate columns keeps activities finite", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("flashier.utils")
  set.seed(9)
  n <- 30; g <- 30
  Y <- .synth_counts(n = n, g = g, seed = 9)$Y
  meta <- data.frame(observation_id = paste0("s", 1:n),
                     age = seq(20, 60, length.out = n),
                     age_copy = seq(20, 60, length.out = n))   # aliased with age
  cpd <- as_cell_program_data(list(CT1 = Y), sample_metadata = meta)
  # feature_scale declared on BOTH sides: the projection contract refuses
  # undeclared reference units (query compatibility check).
  fit <- canonicalize_programs(fit_celltype_programs(
    cpd, covariates = c("age", "age_copy"), max_factors = 2, seed = 1,
    feature_scale = "logcounts"))
  beta <- fit$fits$CT1$preproc_replay$residualize_coefficients
  expect_true(anyNA(beta))   # lm.fit pivots the aliased column to NA

  nq <- 20
  expr <- Y[1:nq, , drop = FALSE]
  rownames(expr) <- paste0("c", seq_len(nq))
  cell_types <- setNames(rep("CT1", nq), paste0("c", seq_len(nq)))
  donors <- setNames(paste0("s", seq_len(nq)), paste0("c", seq_len(nq)))
  pr <- project_programs_to_cells(fit, expr, cell_types, donors = donors,
                                  sample_metadata = meta,
                                  feature_scale = "logcounts")
  expect_false(anyNA(pr$activity))
  expect_false(any(is.infinite(pr$activity)))
})

# ---- cebmf_torch backend (needs a python env) -------------------------------

test_that("cebmf_torch backend rejects fixed-mode covariates", {
  set.seed(6)
  n <- 30; g <- 25
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  meta <- data.frame(observation_id = paste0("s", 1:n), age = seq(20, 60, length.out = n))
  cpd <- as_cell_program_data(list(CT1 = Y), sample_metadata = meta)
  expect_error(fit_celltype_programs(cpd, backend = "cebmf_torch",
                                     covariates = ~age, covariate_mode = "fixed"),
               "no fixed-loadings path", fixed = TRUE)
})

test_that("cebmf_torch backend fits, canonicalizes, and reports NA elbo", {
  skip_if(!cebmf_env_ready(), "no python env for cebmf_torch backend")
  Y <- .synth_counts()$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  fit <- fit_celltype_programs(cpd, backend = "cebmf_torch", max_factors = 4, seed = 1)
  fl <- fit$fits$CT1$flash
  expect_equal(fl$backend, "cebmf_torch")
  expect_true(fl$n_factors >= 1)
  expect_true(is.matrix(fl$L_pm) && nrow(fl$L_pm) == nrow(Y))
  expect_equal(ncol(fl$L_pm), fl$n_factors)
  # registered mapping contract: nothing manufactured
  expect_null(fl$elbo)
  expect_null(fl$L_psd); expect_null(fl$F_psd); expect_null(fl$F_lfsr)
  expect_true(all(is.finite(fl$objective_trace)))
  fit <- canonicalize_programs(fit)
  expect_true(ncol(fit$scores$CT1) >= 1)
  summ <- program_summary(fit)
  expect_true(is.na(summ$elbo[summ$cell_type == "CT1"]))
})

test_that("cebmf_torch backend is deterministic at a fixed seed", {
  skip_if(!cebmf_env_ready(), "no python env for cebmf_torch backend")
  Y <- .synth_counts()$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  f1 <- fit_celltype_programs(cpd, backend = "cebmf_torch", max_factors = 4, seed = 42)
  f2 <- fit_celltype_programs(cpd, backend = "cebmf_torch", max_factors = 4, seed = 42)
  expect_true(all.equal(f1$fits$CT1$flash$L_pm,
                        f2$fits$CT1$flash$L_pm) == TRUE)
  expect_true(all.equal(f1$fits$CT1$flash$F_pm,
                        f2$fits$CT1$flash$F_pm) == TRUE)
})

test_that("cebmf_torch backend recovers the planted programs (parity bar 0.8)", {
  skip_if(!cebmf_env_ready(), "no python env for cebmf_torch backend")
  skip_if_not_installed("flashier")
  Y <- .synth_counts()$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  fc <- fit_celltype_programs(cpd, backend = "cebmf_torch", max_factors = 4, seed = 1)
  ff <- fit_celltype_programs(cpd, backend = "flashier", max_factors = 4, seed = 1)
  Lc <- fc$fits$CT1$flash$L_pm
  Lf <- ff$fits$CT1$flash$L_pm
  k <- min(ncol(Lc), ncol(Lf))
  cors <- abs(cor(Lc, Lf))
  # greedy max-matched absolute correlation
  matched <- numeric(k)
  for (j in seq_len(k)) {
    idx <- which(cors == max(cors), arr.ind = TRUE)[1, ]
    matched[j] <- cors[idx[1], idx[2]]
    cors[idx[1], ] <- -1; cors[, idx[2]] <- -1
  }
  expect_gte(min(matched), 0.8)
})

test_that("var_floor passes through the cebmf_torch backend as a fixed S", {
  skip_if(!cebmf_env_ready(), "no python env for cebmf_torch backend")
  Y <- .synth_counts()$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  fit <- fit_celltype_programs(cpd, backend = "cebmf_torch", max_factors = 4,
                               seed = 1, var_floor = 0.5)
  expect_identical(fit$fits$CT1$S_used, 0.5)
  expect_identical(fit$preprocessing$var_floor, 0.5)
  expect_true(fit$fits$CT1$flash$n_factors >= 1)
})

# ---- var_floor on the flashier backend --------------------------------------

test_that("var_floor numeric scalar fixes S and is recorded in provenance", {
  skip_if_not_installed("flashier")
  Y <- .synth_counts(n = 30, g = 40)$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  fit <- fit_celltype_programs(cpd, max_factors = 3, seed = 1, var_floor = 0.5)
  expect_identical(fit$fits$CT1$S_used, 0.5)
  expect_identical(fit$preprocessing$var_floor, 0.5)
  expect_identical(fit$provenance$var_floor, 0.5)
})

test_that("var_floor vignette computes the documented deterministic S", {
  skip_if_not_installed("flashier")
  Y <- .synth_counts(n = 30, g = 40)$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  fit <- fit_celltype_programs(cpd, max_factors = 3, seed = 1, var_floor = "vignette")
  n <- nrow(fit$fits$CT1$Y_used)
  set.seed(20261016)   # the registered fixed local seed
  s_expected <- sd(log(stats::rpois(1e7, 1 / n) + 1))
  expect_identical(fit$fits$CT1$S_used, s_expected)
  expect_identical(fit$provenance$var_floor, "vignette")
})

test_that("var_floor reroutes the flashier.utils fast path through legacy flash with a warning", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("flashier.utils")
  Y <- .synth_counts(n = 30, g = 40)$Y
  cpd <- as_cell_program_data(list(CT1 = Y))
  fit <- NULL
  expect_warning(fit <- fit_celltype_programs(cpd, max_factors = 3, seed = 1,
                                              var_floor = 0.5),
                 "legacy flashier backend")
  expect_identical(fit$fits$CT1$S_used, 0.5)
})

# ---- edge cases: basilisk env spec + python option handling -----------------

test_that("cebmf basilisk env spec satisfies the pip invariants CI taught us", {
  skip_if_not(requireNamespace("basilisk", quietly = TRUE))
  env <- .cebmf_basilisk_env()
  # basilisk <= 1.2x keeps a separate `pip` slot; newer basilisk deprecates the
  # constructor arg and APPENDS pip entries to `packages` (slot removed). The
  # conda spec here is exactly c("python=3.11"), so setdiff recovers the pip
  # entries under both layouts.
  pip <- if ("pip" %in% slotNames(env)) env@pip else setdiff(env@packages, "python=3.11")
  # basilisk 1.18 passes pip entries through system2 unquoted (shell splits
  # on whitespace); basilisk 1.24 shell-quotes each entry (a flag entry
  # becomes a quoted "package name"). Both break on flag entries or spaces.
  expect_false(any(startsWith(pip, "-")),
               info = paste("flag entries break basilisk 1.24:", pip[startsWith(pip, "-")]))
  expect_false(any(grepl("\\s", pip)),
               info = paste("whitespace breaks basilisk 1.18 system2:", pip[grepl("\\s", pip)]))
  # every entry is either ==-pinned or a PEP 508 direct reference
  ok <- grepl("==", pip) | grepl("^[A-Za-z0-9_.-]+@git\\+", pip)
  expect_true(all(ok), info = paste("unpinned/unparseable:", pip[!ok]))
  # python itself is pinned in the conda/packages spec
  expect_true(any(grepl("^python=", env@packages)))
  # the direct reference stays commit-pinned (40 hex chars)
  git_ref <- grep("@git\\+", pip, value = TRUE)
  expect_true(all(grepl("@[0-9a-f]{40}$", git_ref)),
              info = "cebmf_torch must be pinned to a full commit SHA")
})

test_that(".cebmf_use_python_env aborts cleanly on bad paths", {
  skip_if_not(requireNamespace("reticulate", quietly = TRUE))
  expect_error(.cebmf_use_python_env("/definitely/not/here"),
               "non-existent python", fixed = TRUE)
  # an existing directory WITHOUT bin/python must abort here too, not fall
  # through to reticulate's use_python error (file.exists is TRUE for dirs)
  d <- tempfile(); dir.create(d)
  expect_error(.cebmf_use_python_env(d), "non-existent python", fixed = TRUE)
})
