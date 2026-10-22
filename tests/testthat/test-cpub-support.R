# CPUB-009/010/012/027 bounded support contracts. Also runnable by sourcing
# package R files and using testthat::test_file(..., env = that_environment).

.support_block_fit <- function() {
  ids <- paste0("d", 1:8)
  make <- function(ct) matrix(seq_len(16), 8, 2,
                             dimnames = list(ids, paste0(ct, c("g1", "g2"))))
  list(cell_types = c("A", "B"),
       fits = list(A = list(Y_used = make("A"), flash = list(n_factors = 1L)),
                   B = list(Y_used = make("B"), flash = list(n_factors = 1L))),
       scores = list(A = matrix(seq_len(8), 8, 1, dimnames = list(ids, "A_1")),
                     B = matrix(seq_len(8), 8, 1, dimnames = list(ids, "B_1"))),
       loadings = list(A = matrix(c(1, 2), 2, 1), B = matrix(c(1, 2), 2, 1)))
}

test_that("block level/order x K0/K1 uses explicit support and honest membership", {
  for (k in c(0L, 1L)) {
    calls <- 0L
    testthat::local_mocked_bindings(
      flash = function(data, ...) {
        calls <<- calls + 1L
        list(n_factors = k,
             L_pm = matrix(1, nrow(data), k),
             F_pm = matrix(1 / sqrt(ncol(data)), ncol(data), k))
      }, .package = "flashier")
    for (level in c("scores", "expression")) for (order in c("shared_first", "private_first")) {
      fit <- .support_block_fit()
      before <- calls
      if (level == "scores" && order == "private_first") {
        # A global R must never be used as an accidental solver input.
        e <- new.env(parent = environment(fit_joint_block))
        e$R <- matrix(999, 8, 2)
        f <- fit_joint_block
        environment(f) <- e
        expect_error(f(fit, level, order), "unsupported")
        expect_identical(calls, before)
        next
      }
      b <- fit_joint_block(fit, level, order)
      expect_equal(dim(b$block_mass), c(k, 2L))
      expect_equal(colnames(b$block_mass), c("A", "B"))
      cl <- block_clusters(b)
      if (k) {
        expect_equal(as.numeric(b$block_mass), c(0.5, 0.5))
        expect_setequal(cl$clusters$cell_type, c("A", "B"))
        expect_true(all(is.na(cl$clusters$fwer_p)))
        expect_equal(rownames(cl$activity), paste0("d", 1:8))
        if (level == "expression") {
          expect_identical(rownames(b$gene_loadings$A), c("Ag1", "Ag2"))
        }
      } else {
        expect_equal(nrow(cl$clusters), 0L)
        expect_null(cl$activity)
      }
    }
  }
})

test_that("old scores blocks without membership fail rather than silently extract empty calls", {
  b <- list(shared_fit = list(n_factors = 1L),
            summary = data.frame(class = "partial:A+B", factor_id = "SF01"))
  expect_error(block_clusters(b), "block_mass")
})

test_that("observed-stratum QR equals a hand-centered rank oracle", {
  ids <- paste0("d", 1:12)
  blocks <- setNames(rep(c("a", "b", "c"), each = 4), ids)
  S <- cbind(p1 = seq_len(12), p2 = c(NA_real_, NA_real_, 3:12),
             p3 = c(rep(NA_real_, 8), 1:4))
  rownames(S) <- ids
  S[9:12, "p1"] <- NA_real_ # entire stratum c absent for p1
  z <- .cmp_prep(list(A = S), blocks, min_donors = 4)
  r <- .cmp_rank_gauss(S)
  for (j in 1:3) {
    obs <- !is.na(S[, j])
    oracle <- r[obs, j] - ave(r[obs, j], blocks[ids[obs]])
    expect_equal(as.numeric(z$A[obs, j]), as.numeric(oracle), tolerance = 1e-12)
    expect_identical(is.na(z$A[, j]), is.na(S[, j]))
  }
  info <- attr(z, "program_support")
  expect_equal(info$design_rank, c(2L, 3L, 1L))
  expect_equal(info$residual_df, c(6L, 7L, 3L))
  expect_true(all(info$callable))
})

test_that("complete full-rank projection stays equivalent and invalid residuals are reported", {
  ids <- paste0("d", 1:10)
  blocks <- setNames(rep(c("a", "b"), each = 5), ids)
  S <- cbind(varying = c(3, 1, 5, 4, 2, 8, 7, 10, 6, 9),
             constant = rep(1, 10), missing = rep(NA_real_, 10),
             sparse = c(1, 2, rep(NA_real_, 8)))
  rownames(S) <- ids
  z <- .cmp_prep(list(A = S), blocks, min_donors = 4)
  r <- .cmp_rank_gauss(S)[, 1]
  X <- model.matrix(~ factor(blocks))
  expect_equal(as.numeric(z$A[, 1]), as.numeric(r - X %*% qr.solve(X, r)), tolerance = 1e-12)
  info <- attr(z, "program_support")
  expect_equal(info$reason, c("ok", "constant_residual", "all_missing", "insufficient_donors"))
  expect_true(all(is.na(z$A[, "missing"])))
  other <- S[, 1, drop = FALSE]; colnames(other) <- "other"
  zz <- .cmp_prep(list(A = S, B = other), blocks, min_donors = 4)
  ee <- .cmp_edges(zz, .cmp_pairs(zz, 4))
  expect_false(ee$callable[ee$prog_a == "missing"])
  expect_true(is.na(ee$stat[ee$prog_a == "missing"]))
  missing_ct <- .cmp_prep(list(A = S[, "missing", drop = FALSE]), blocks, min_donors = 4)
  expect_length(missing_ct, 0L)
  expect_equal(attr(missing_ct, "program_support")$reason, "all_missing")
  sat <- .cmp_prep(list(A = S[, 1, drop = FALSE]), setNames(ids, ids), min_donors = 4)
  expect_equal(attr(sat, "program_support")$reason, "no_residual_df")
})

test_that("overlap thresholds and disconnected support never become scientific zero edges", {
  ids <- paste0("d", 1:12)
  A <- matrix(c(1:6, rep(NA_real_, 6)), 12, 1, dimnames = list(ids, "A_1"))
  B <- matrix(c(rep(NA_real_, 6), 1:6), 12, 1, dimnames = list(ids, "B_1"))
  Z <- .cmp_prep(list(A = A, B = B), setNames(rep("all", 12), ids), min_donors = 4)
  E <- .cmp_edges(Z, .cmp_pairs(Z, 5))
  expect_equal(E$n, 0)
  expect_false(E$callable)
  expect_true(is.na(E$stat))
  expect_equal(E$reason, "insufficient_shared")
  B <- A; colnames(B) <- "B_1"; B[5:6, ] <- NA_real_
  Z <- .cmp_prep(list(A = A, B = B), setNames(rep("all", 12), ids), min_donors = 4)
  E <- .cmp_edges(Z, .cmp_pairs(Z, 5))
  expect_equal(E$n, 4)
  expect_false(E$callable)
  expect_true(is.na(E$stat))
})

test_that("all-constant CMP nodes expose NA p values and support receipts", {
  ids <- paste0("d", 1:8)
  M <- matrix(1, 8, 1, dimnames = list(ids, "A_1"))
  N <- M; colnames(N) <- "B_1"
  fit <- list(scores = list(A = M, B = N), loadings = list(A = matrix(1), B = matrix(1)),
              observation_ids = ids)
  res <- find_coordinated_programs(fit, min_donors = 4, min_shared = 4,
                                    n_perm = 20, ambient = "none")
  expect_equal(nrow(res$cmps), 0L)
  expect_true(all(is.na(res$nodes$node_p)))
  expect_true(all(is.na(res$nodes$node_stat)))
  expect_false(any(res$nodes$callable))
  expect_false(any(res$edge_support$callable))
  expect_true(all(res$program_support$reason == "constant_residual"))
})

.support_stability_fit <- function(backend = "flashier", floor = "none", spec = NULL) {
  ids <- paste0("d", 1:8)
  Y <- matrix(seq_len(16), 8, 2, dimnames = list(ids, c("g1", "g2")))
  list(pseudobulk = list(A = Y), cell_types = "A", shared_genes = colnames(Y),
       observation_ids = ids, sample_metadata = data.frame(observation_id = ids),
       fits = list(A = list()), canonicalization = list(),
       loadings = list(A = matrix(1, 2, 1, dimnames = list(colnames(Y), "A_1"))),
       scores = list(A = matrix(seq_len(8), 8, 1, dimnames = list(ids, "A_1"))),
       preprocessing = list(loading_prior = "point_normal", max_factors = 1L,
                            var_type = 2L, backfit = FALSE, nullcheck = FALSE,
                            centered = TRUE, scaled = FALSE,
                            covariates = ~ age, covariate_mode = "residualize",
                            flash_control = list(greedy_Kmax = 1L),
                            feature_scale = "log_expression",
                            feature_sets = list(A = c("g1", "g2")),
                            backend = backend, var_floor = floor, spec = spec),
       provenance = list(seed = 17L))
}

test_that("stability captures estimator controls on frozen and legacy paths including failures", {
  for (frozen in c(FALSE, TRUE)) for (backend in c("flashier", "cebmf_torch")) {
    for (floor in list("none", "vignette", 0.2)) {
      spec <- if (frozen) list(digest = "fixture", feature_scale = "log_expression") else NULL
      fit <- .support_stability_fit(backend, floor, spec)
      seen <- list()
      e <- new.env(parent = environment(assess_program_stability))
      e$fit_celltype_programs <- function(...) {
        args <- list(...); seen[[length(seen) + 1L]] <<- args
        if (length(seen) == 2L) stop("intentional refit failure")
        list(scores = list(A = matrix(numeric(0), nrow(args[[1]]$pseudobulk$A), 0)))
      }
      e$canonicalize_programs <- identity
      assess <- assess_program_stability; environment(assess) <- e
      expect_warning(out <- assess(fit, n_boot = 3, progress = FALSE, seed = 3),
                     "1 of 3 subsample refits failed")
      expect_length(seen, 3L)
      expect_equal(out$stability$valid_refits, 2L)
      expect_equal(out$stability$failed_refits, 1L)
      expect_equal(out$stability$replicate_log$status, c("ok", "error", "ok"))
      expect_match(out$stability$replicate_log$error_message[2], "intentional")
      expect_equal(out$stability$per_program$n_eligible, 0L)
      for (b in 1:3) {
        args <- seen[[b]]
        expect_identical(args$backend, backend)
        expect_identical(args$var_floor, floor)
        expect_identical(args$loading_prior, "point_normal")
        expect_identical(args$flash_control, fit$preprocessing$flash_control)
        expect_identical(out$stability$replicate_log$refit_settings[[b]], args[-1L])
        if (frozen) expect_identical(args$preprocessing_spec, spec) else {
          expect_identical(args$feature_scale, "log_expression")
          expect_identical(args$features, fit$preprocessing$feature_sets)
          expect_identical(args$covariates, fit$preprocessing$covariates)
        }
      }
    }
  }
})

test_that("stability all-failure denominator is NA and legacy defaults are explicit", {
  fit <- .support_stability_fit()
  fit$preprocessing$backend <- fit$preprocessing$var_floor <- NULL
  a <- .stability_refit_args(list(), fit$preprocessing, 1L)
  expect_identical(a$backend, "flashier")
  expect_identical(a$var_floor, "none")
  e <- new.env(parent = environment(assess_program_stability))
  e$fit_celltype_programs <- function(...) stop("always fails")
  assess <- assess_program_stability; environment(assess) <- e
  expect_warning(out <- assess(fit, n_boot = 2, progress = FALSE), "all 2")
  expect_true(is.na(out$stability$per_program$recovery_freq))
  expect_true(is.na(out$stability$per_program$conditional_recovery_freq))
  expect_equal(nrow(out$stability$replicate_log), 2L)
})

test_that("tiny real flashier stability refits preserve numeric fixed-noise settings", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  old <- options(cellprograms.use_flashier_utils = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(91)
  ids <- paste0("d", 1:24)
  Y <- outer(rnorm(24), c(rep(2, 8), rep(-2, 8))) + matrix(rnorm(24 * 16, sd = 0.1), 24)
  dimnames(Y) <- list(ids, paste0("g", 1:16))
  x <- as_cell_program_data(list(A = Y))
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 1L,
                              var_floor = 0.2, seed = 1L))
  expect_identical(fit$preprocessing$var_floor, 0.2)
  out <- assess_program_stability(fit, n_boot = 1L, seed = 2L, progress = FALSE)
  expect_equal(out$stability$valid_refits, 1L)
  expect_identical(out$stability$replicate_log$refit_settings[[1]]$var_floor, 0.2)
  expect_identical(out$stability$replicate_log$refit_settings[[1]]$backend, "flashier")
  eff <- out$stability$replicate_log$effective_refit[[1]]
  expect_true(all(vapply(eff, function(f) f$effective_config$noise_mode ==
    "known_fixed" && f$S_used == .2, logical(1))))
  expect_identical(out$stability$target,
    "full_fitting_procedure_including_documented_fallback_ladder")
})

.support_seurat <- function() {
  set.seed(7)
  counts <- matrix(rpois(80 * 32, 3), 80, 32,
                    dimnames = list(paste0("g", 1:80), paste0("c", 1:32)))
  Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
}

test_that("automatic PCA returns changed object and consistent lowercase reduction key", {
  skip_if_not_installed("Seurat")
  obj <- .support_seurat()
  red <- .sc_prepare_reduction(obj, NULL, 1:2, 40, 3, 7)
  expect_false("pca" %in% names(obj@reductions))
  expect_true("pca" %in% names(red$object@reductions))
  expect_identical(red$reduction, "pca")
  expect_identical(red$dims, 1:2)
  expect_equal(rownames(red$object@reductions$pca@cell.embeddings), colnames(obj))
  for (key in c("pca", "PCA")) {
    again <- .sc_prepare_reduction(red$object, key, 1:2, 40, 3, 7)
    expect_identical(again$reduction, "pca")
    expect_identical(again$object, red$object)
  }
  expect_error(.sc_prepare_reduction(red$object, "pca", 10, 40, 3, 7), "dims")
  expect_error(.sc_prepare_reduction(red$object, "missing", 1, 40, 3, 7), "not found")
})

test_that("automatic reduction grouped SuperCell integration conserves count sums", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SuperCell")
  obj <- .support_seurat()
  obj$donor <- rep(c("d1", "d2"), each = 16)
  out <- MetacellsByGroupsSuperCell(obj, group.by = "donor", reduction = NULL,
                                    n_pc = 3, n_top_genes = 40, min_cells = 16,
                                    k.knn = 5, cells_per_metacell = 4,
                                    run_mcrigor = FALSE, verbose = FALSE)
  mc <- GetMetacellSuperCell(out)$metacell_object
  expect_equal(sum(SeuratObject::LayerData(mc, layer = "counts")),
               sum(SeuratObject::LayerData(obj, layer = "counts")))
  expect_true(all(grepl("^donor_", colnames(mc))))
  expect_setequal(mc$donor, c("d1", "d2"))
})
