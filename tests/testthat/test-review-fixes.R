# Fixes from the external scientific-correctness review (Biomni, 2026-10-08)

test_that("fixed-mode design survives a factor with an unused level", {
  meta <- data.frame(observation_id = paste0("s", 1:20),
                     site = factor(rep(c("a", "b"), 10), levels = c("a", "b", "c")))
  X <- .cov_design(meta$observation_id, meta, "site", "CT", intercept = FALSE)
  expect_equal(ncol(X), 1L)
  expect_false(anyNA(X))
})

test_that("a level with a single donor triggers a leverage warning; too few residual df aborts", {
  meta <- data.frame(observation_id = paste0("s", 1:30),
                     g = c(rep("a", 14), rep("b", 15), "c"))
  expect_warning(.cov_design(meta$observation_id, meta, "g", "CT"), "leverage")
  small <- data.frame(observation_id = paste0("s", 1:12), g = rep(letters[1:6], 2))
  expect_error(.cov_design(small$observation_id, small, "g", "CT"), "residual degrees of freedom")
})

test_that("per-pair seeds differ between pairs and are reproducible", {
  expect_false(.pair_seed(1, c("A", "B")) == .pair_seed(1, c("A", "C")))
  expect_false(.pair_seed(1, c("A", "B")) == .pair_seed(1, c("B", "A")))
  expect_identical(.pair_seed(1, c("A", "B")), .pair_seed(1, c("A", "B")))
})

test_that("sharing_spectrum scales n_perm with the number of pairs and warns when too small", {
  set.seed(2); n <- 30
  mk <- function(ct) matrix(rnorm(n * 2), n, 2, dimnames = list(paste0("s", 1:n), paste0(ct, "_", 1:2)))
  sc <- list(A = mk("A"), B = mk("B"), C = mk("C"))
  fit <- structure(list(scores = sc, loadings = lapply(sc, function(S) S[0, ])), class = "cell_program_fit")
  expect_no_warning(r <- sharing_spectrum(fit, seed = 1))
  expect_gte(r$pairs[[1]]$n_perm, 999L)
  expect_warning(sharing_spectrum(fit, n_perm = 50, seed = 1), "permutation floor")
})

test_that("gene weights are capped at 10x the median", {
  set.seed(3)
  Y <- matrix(rnorm(40 * 6), 40, 6, dimnames = list(paste0("s", 1:40), paste0("g", 1:6)))
  Z <- matrix(rnorm(40), 40, 1, dimnames = list(rownames(Y), "A_1"))
  W <- matrix(c(1, rep(0, 5)), 6, 1, dimnames = list(colnames(Y), "A_1"))
  Y[, 1] <- Z[, 1]                         # gene 1 is fitted perfectly -> residual sd ~ 0
  fit <- list(fits = list(A = list(Y_used = Y), B = list(Y_used = Y)),
              scores = list(A = Z, B = Z), loadings = list(A = W, B = W))
  w <- .gene_weights(fit, "A", "B", colnames(Y))
  expect_lte(max(w) / median(w), 10 + 1e-8)
})

test_that("a silent NaN ELBO triggers the fallback ladder", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  real <- flashier::flash
  testthat::local_mocked_bindings(
    flash = function(...) { f <- real(...); if (isTRUE(list(...)$nullcheck)) f$elbo <- NaN; f },
    .package = "flashier")
  fit <- suppressWarnings(fit_celltype_programs(.toy(), max_factors = 2,
                                                loading_prior = "point_normal", seed = 1))
  expect_identical(fit$fits$A$fallback, "no_nullcheck")
})

test_that("degenerate permutation null gives NA z instead of NaN/Inf", {
  one <- matrix(1, 20, 1)   # constant column: every row permutation gives the same statistic, sd(null) = 0
  r <- .pa_engine(one, one, n_perm = 20, seed = 1)
  expect_true(is.na(r$z_pair))
  expect_true(is.finite(r$p_pair))
})

# ---- 2026-10-19 implementation review: criticals ----

test_that("multi-group soft assignments combine block-diagonally (C-1)", {
  mk <- function(cells, mcs, w) {
    Matrix::Matrix(w, length(cells), length(mcs), sparse = TRUE,
                   dimnames = list(cells, mcs))
  }
  a <- mk(paste0("c", 1:3), paste0("grpA_", 1:2), 0.5)
  b <- mk(paste0("c", 4:6), paste0("grpB_", 1:4), 0.25)   # unequal metacell counts
  s <- .combine_soft_assignments(list(A = a, B = b))
  expect_s4_class(s, "dgCMatrix")
  expect_equal(dim(s), c(6L, 6L))
  expect_equal(rownames(s), paste0("c", 1:6))
  expect_equal(colnames(s), c(paste0("grpA_", 1:2), paste0("grpB_", 1:4)))
  # group-local weights: no cross-group leakage
  expect_true(all(s[1:3, 3:6] == 0))
  expect_true(all(s[4:6, 1:2] == 0))
  expect_equal(Matrix::nnzero(s), Matrix::nnzero(a) + Matrix::nnzero(b))
})

test_that("cmp_context_association never reports x* covariates as context terms (C-2)", {
  set.seed(11)
  ids <- paste0("d", seq_len(40))
  z <- rnorm(40)
  mk <- function(ct) {
    S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
    S[, 1] <- z + rnorm(40, sd = 0.3)
    S
  }
  sc <- list(A = mk("A"), B = mk("B"))
  ld <- lapply(sc, function(S) S[0, ])
  meta <- data.frame(observation_id = ids,
                     group = rep(c("ctrl", "stim"), each = 20),
                     xsite = rep(c("p", "q"), 20))
  fit <- structure(list(scores = sc, loadings = ld, cell_types = c("A", "B"),
                        observation_ids = ids, sample_metadata = meta),
                   class = c("cell_program_fit", "cell_program_data", "list"))
  cmp <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
  res <- cmp_context_association(fit, cmp, "group", covariates = "xsite")
  expect_gt(nrow(res), 0)
  expect_true(all(grepl("^\\.ctx", res$term)))
  expect_false(any(grepl("xsite", res$term)))
  expect_error(cmp_context_association(fit, cmp, "group", covariates = "y"), "reserved")
})

# ---- 2026-10-19 implementation review: P1 crash guards ----

test_that("scale=TRUE with a constant gene drops it instead of crashing (M-1)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(6)
  n <- 24; g <- 30
  Y <- matrix(rnorm(n * g), n, g, dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  Y[, "g5"] <- 3.14   # constant gene: NaN under center+scale
  cpd <- as_cell_program_data(list(CT1 = Y))
  expect_message(
    fit <- fit_celltype_programs(cpd, max_factors = 2, center = TRUE, scale = TRUE,
                                 loading_prior = "point_normal", seed = 1),
    "zero-variance")
  expect_false("g5" %in% fit$fits$CT1$genes)
})

test_that("stability resampling warns instead of crashing for <5-donor cell types (M-5)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(7)
  mk <- function(n, g = 25) matrix(rnorm(n * g), n, g,
                                   dimnames = list(paste0("s", 1:n), paste0("g", 1:g)))
  cpd <- as_cell_program_data(list(A = mk(10), B = mk(4)))
  fit <- canonicalize_programs(fit_celltype_programs(cpd, max_factors = 2,
                                                     loading_prior = "point_normal", seed = 1))
  expect_warning(
    fit <- assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE),
    "fewer than 5 donors")
  expect_true(!is.null(fit$stability))
})

test_that("plot_program_scores handles per-cell-type observation sets (M-6)", {
  set.seed(9)
  sa <- matrix(rnorm(20), 10, 2, dimnames = list(paste0("s", 1:10), c("A_1", "A_2")))
  sb <- matrix(rnorm(6), 6, 1, dimnames = list(paste0("s", 5:10), "B_1"))
  fit <- structure(list(scores = list(A = sa, B = sb),
                        sample_metadata = data.frame(observation_id = paste0("s", 1:10))),
                   class = "cell_program_fit")
  expect_no_error(p <- plot_program_scores(fit))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_sharing_spectrum handles >9 pairs and labels shared pairs with q (M-7)", {
  set.seed(10)
  pairs <- t(combn(LETTERS[1:5], 2))   # 10 pairs > brewer Set1's 9-color limit
  s <- data.frame(min_angle_deg = runif(10, 0, 90),
                  pair = paste(pairs[, 1], pairs[, 2], sep = "-"),
                  underpowered = FALSE,
                  shared = rep(c(TRUE, FALSE), 5),
                  q_pair = runif(10, 0, 0.2),
                  n_shared_05 = sample(0:3, 10, TRUE))
  p <- plot_sharing_spectrum(list(summary = s))
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_true(any(grepl("q = ", levels(p$data$pair_lab))))
})

test_that("plot_artifact_flags works for any screened set and CMP-shaped input (M-8)", {
  fl <- data.frame(cell_type = c("A", "B"), program_id = c("A_1", "B_1"),
                   sex_frac = c(0.5, 0), sex_flag = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)
  p <- plot_artifact_flags(fl)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_true("sex" %in% p$data$set)
  flc <- data.frame(cmp = 1L, n_members = 2L, vdj_frac = 0.25, vdj_flag = TRUE)
  expect_s3_class(plot_artifact_flags(flc), "ggplot")
})

.union_toy <- function() {
  set.seed(8)
  ids <- paste0("s", 1:12)
  Y <- matrix(rnorm(12 * 6), 12, 6, dimnames = list(ids, paste0("g", 1:6)))
  Y[, "g3"] <- 7   # constant gene: NaN under center+scale
  S <- matrix(rnorm(12 * 2), 12, 2, dimnames = list(ids, c("A_1", "A_2")))
  structure(list(
    fits = list(A = list(genes = paste0("g", 1:6), flash = list(fixed_L = NULL))),
    scores = list(A = S),
    pseudobulk = list(A = Y),
    preprocessing = list(covariates = NULL, covariate_mode = "residualize",
                         centered = TRUE, scaled = TRUE),
    sample_metadata = data.frame(observation_id = ids)),
    class = "cell_program_fit")
}

test_that("union_panel_loadings leaves constant panel genes NA instead of aborting (M-9)", {
  fit <- .union_toy()
  up <- union_panel_loadings(fit)
  expect_true(all(is.na(up$loadings$A["g3", ])))
  expect_true(all(is.finite(up$loadings$A[setdiff(paste0("g", 1:6), "g3"), ])))
  # zero-overlap panel genes stay NA, no error
  up0 <- union_panel_loadings(fit, genes = c("g1", "nope"))
  expect_true(all(is.na(up0$loadings$A["nope", ])))
  expect_true(all(is.finite(up0$loadings$A["g1", ])))
  up1 <- union_panel_loadings(fit, genes = "nope")
  expect_true(all(is.na(up1$loadings$A)))
})

test_that("union_panel_loadings replays a rank-deficient fixed design without solve() crash", {
  fit <- .union_toy()
  fit$preprocessing$covariates <- "batch"
  fit$preprocessing$covariate_mode <- "fixed"
  Xc <- cbind(rep(1, 12), rep(c(0, 1), 6), rep(c(0, 1), 6))  # col2 == col3: singular crossprod
  fit$fits$A$flash$fixed_L <- Xc
  up <- union_panel_loadings(fit)
  expect_true(all(is.finite(up$loadings$A["g1", ])))
})

# ---- 2026-10-19 implementation review: P2 quality/metadata ----

test_that("cebmf_torch backend warns on non-converged fits (M-2)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  set.seed(12)
  Y <- matrix(rnorm(24 * 30), 24, 30, dimnames = list(paste0("s", 1:24), paste0("g", 1:30)))
  cpd <- as_cell_program_data(list(CT1 = Y))
  fake <- list(
    n_factors = 2L, elbo = NULL, objective_trace = c(-100, -90),
    objective_direction = "maximize", residuals_sd = NULL,
    L_pm = matrix(rnorm(24 * 2), 24, 2), F_pm = matrix(rnorm(30 * 2), 30, 2),
    L_psd = NULL, F_psd = NULL, F_lfsr = NULL, fixed_L = NULL, fixed_F = NULL,
    tau = 1, effective_k = 2L, converged = FALSE, guard_diagnostics = list(),
    backend = "cebmf_torch"
  )
  testthat::local_mocked_bindings(.fit_via_cebmf = function(...) fake)
  expect_warning(
    fit_celltype_programs(cpd, max_factors = 2, backend = "cebmf_torch", seed = 1),
    "did not converge")
})

test_that("project_programs_to_cells refuses a mismatched feature_scale (M-4)", {
  skip_if_not_installed("flashier.utils")
  fit <- structure(list(
    scores = list(A = matrix(0, 2, 1, dimnames = list(c("s1", "s2"), "A_1"))),
    canonicalization = list(),
    fits = list(A = list(flash = structure(list(), flashbridge_metadata = list()),
                         preproc_replay = list(center = 0, scale = 1,
                                               residualize_coefficients = NULL,
                                               design_columns = NULL),
                         genes = "g1")),
    preprocessing = list(feature_scale = "logcounts", covariates = NULL),
    cell_types = "A"),
    class = "cell_program_fit")
  expr <- matrix(0, 2, 3, dimnames = list(c("c1", "c2"), c("g1", "g2", "g3")))
  cts <- c(c1 = "A", c2 = "A")
  expect_error(
    project_programs_to_cells(fit, expr, cts, feature_scale = "counts"),
    "differs from the fit's recorded units")
  # a matching declaration passes the units gate (any later error is fine)
  err <- tryCatch(
    project_programs_to_cells(fit, expr, cts, feature_scale = "logcounts"),
    error = function(e) e)
  expect_true(inherits(err, "error"))
  expect_false(grepl("differs from the fit's recorded units", conditionMessage(err)))
})

# ---- 2026-10-19 implementation review: P3 batch 3 (plots / match_programs / mofacellular) ----

test_that("plot_context_association omits non-finite CIs and labels .ctx plainly (43)", {
  assoc <- data.frame(cmp = c(1L, 1L), term = c(".ctx", "batch"),
                      estimate = c(0, 0.5), t = c(0, 2.5),
                      p = c(1, 0.05), n = c(20L, 20L))
  expect_warning(p <- plot_context_association(assoc), "zero")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_true("context variable" %in% p$data$term)
  assoc_f <- data.frame(cmp = 1L, term = ".ctxhigh", estimate = 0.4, t = 2,
                        p = 0.06, n = 20L)
  pf <- plot_context_association(assoc_f)
  expect_true("context variable: high" %in% pf$data$term)
})

test_that("plot_variance_explained maps PVE via native_indices, never positionally (44)", {
  set.seed(11)
  S <- matrix(rnorm(20), 10, 2, dimnames = list(paste0("s", 1:10), c("A_1", "A_2")))
  # flash object holds 3 factors (intercept at native position 1); the
  # canonicalized scores keep only the 2 free factors.
  fl <- structure(list(pve = c(0.5, 0.3, 0.1)),
                  flashbridge_metadata = list(native_indices = 2:3))
  fit <- structure(list(fits = list(A = list(flash = fl)), scores = list(A = S)),
                   class = "cell_program_fit")
  p <- plot_variance_explained(fit)
  expect_equal(p$data$pve, c(0.3, 0.1))
  # length mismatch without a usable mapping -> warn and skip, never mislabel
  fit2 <- structure(list(fits = list(A = list(flash = list(pve = c(0.5, 0.3, 0.1)))),
                         scores = list(A = S)), class = "cell_program_fit")
  expect_warning(p2 <- plot_variance_explained(fit2), "does not match")
  expect_s3_class(p2, "ggplot")
  # uncanonicalized fit aborts loudly instead of an opaque ncol(NULL) error
  fit3 <- structure(list(fits = list(A = list(flash = fl))), class = "cell_program_fit")
  expect_error(plot_variance_explained(fit3), "canonicalize")
})

test_that("plot_program_scores fails closed on unjoinable annotation (45)", {
  set.seed(13)
  S <- matrix(rnorm(20), 10, 2, dimnames = list(paste0("s", 1:10), c("A_1", "A_2")))
  fit <- structure(list(scores = list(A = S),
                        sample_metadata = data.frame(observation_id = paste0("s", 1:10))),
                   class = "cell_program_fit")
  # data.frame annotation without observation_id -> abort (was silent all-NA)
  expect_error(plot_program_scores(fit, annotation = data.frame(sev = 1:2)),
               "observation_id")
  # character annotation naming a missing column -> abort
  expect_error(plot_program_scores(fit, annotation = "nope"), "not found")
  # zero matching observation ids -> warn, still returns a plot
  expect_warning(p <- plot_program_scores(
    fit, annotation = data.frame(observation_id = paste0("x", 1:10), sev = 1:10)),
    "observation_id")
  expect_s3_class(p, "ggplot")
})

test_that("match_programs space = \"scores\" requires shared observations (17)", {
  set.seed(14)
  Wa <- matrix(rnorm(40), 20, 2, dimnames = list(paste0("g", 1:20), c("B_1", "B_2")))
  Wb <- matrix(rnorm(40), 20, 2, dimnames = list(paste0("g", 1:20), c("B_1", "B_2")))
  mk <- function(W, obs) structure(list(
    loadings = list(B = W),
    scores = list(B = matrix(rnorm(length(obs) * ncol(W)), length(obs), ncol(W),
                             dimnames = list(obs, colnames(W))))),
    class = "cell_program_fit")
  fa <- mk(Wa, paste0("s", 1:5)); fb <- mk(Wb, paste0("t", 1:5))
  expect_error(match_programs(fa, fb, "B", space = "scores"), "shared observations")
  fb2 <- mk(Wb, paste0("s", 1:5))
  expect_no_error(match_programs(fa, fb2, "B", space = "scores"))
})

test_that("as_mofacellular_input fails closed on ambiguous names, duplicate keys, empty types (48-50)", {
  g <- paste0("g", 1:5)
  mk <- function(mats) structure(list(pseudobulk = mats, cell_types = names(mats)),
                                 class = "cell_program_data")
  # 48: "B_cell" + "1" collides with "B" + "cell_1" under paste0(ct, "_", donor)
  x <- mk(list(B_cell = matrix(1:5, 1, 5, dimnames = list("1", g)),
               B = matrix(1:5, 1, 5, dimnames = list("cell_1", g))))
  expect_error(as_mofacellular_input(x), "not unique")
  # 49: duplicate cell_counts keys (match() would silently take the first)
  x2 <- mk(list(B = matrix(1:10, 2, 5, dimnames = list(c("d1", "d2"), g))))
  cc <- data.frame(observation_id = c("d1", "d1"), cell_type = c("B", "B"),
                   cell_counts = c(10, 20))
  expect_error(as_mofacellular_input(x2, cc), "duplicate")
  # 50: zero-observation cell type
  x3 <- mk(list(B = matrix(1:10, 2, 5, dimnames = list(c("d1", "d2"), g)),
                T = matrix(numeric(0), 0, 5, dimnames = list(NULL, g))))
  expect_error(as_mofacellular_input(x3), "zero pseudobulk observations")
  # happy path unchanged
  m <- as_mofacellular_input(x2)
  expect_equal(colnames(m$counts), c("B_d1", "B_d2"))
  expect_equal(m$coldata$donor_id, c("d1", "d2"))
})
