# adjust_batches (ROADMAP small-batch option) + span_alignment (tier-2 subspace metrics)

test_that("adjust_batches applies ComBat before covariate adjustment and warns", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  skip_if_not_installed("sva")
  x <- .toy()
  x$sample_metadata$plate <- factor(rep(c("p1", "p2", "p3", "p4"), length.out = 40))
  expect_warning(
    fit <- fit_celltype_programs(x, adjust_batches = "plate", max_factors = 2, seed = 1),
    "residual batch")
  expect_identical(fit$preprocessing$batch_adjusted, "plate")
  # ComBat equalizes batch location: per-batch column means of the adjusted
  # panel are far closer together than in the raw pseudobulk
  bt <- x$sample_metadata$plate[match(rownames(x$pseudobulk$A), x$sample_metadata$observation_id)]
  spread <- function(M) {
    ms <- tapply(seq_len(nrow(M)), bt, function(i) colMeans(M[i, , drop = FALSE]))
    do.call(rbind, ms)
    max(abs(apply(do.call(rbind, ms), 2, function(v) max(v) - min(v))))
  }
  raw <- x$pseudobulk$A[, fit$fits$A$genes, drop = FALSE]
  expect_lt(spread(fit$fits$A$Y_used), 0.1 * spread(raw))
})

test_that("adjust_batches validates its argument and refuses frozen specs", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  x <- .toy()
  expect_error(fit_celltype_programs(x, adjust_batches = "nope", max_factors = 2), "sample_metadata")
  spec <- freeze_preprocessing(x, features = "all")
  expect_error(fit_celltype_programs(x, adjust_batches = "batch", max_factors = 2,
                                     preprocessing_spec = spec),
               "frozen")
})

test_that("projection refuses batch-adjusted fits (fail-closed)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  skip_if_not_installed("sva")
  x <- .toy()
  fit <- canonicalize_programs(fit_celltype_programs(x, adjust_batches = "batch",
                                                     max_factors = 2, seed = 1,
                                                     feature_scale = "logcounts"))
  ids_a <- rownames(x$pseudobulk$A)
  expect_error(
    project_programs_to_cells(fit, x$pseudobulk$A,
                              cell_types = setNames(rep("A", length(ids_a)), ids_a),
                              donors = setNames(ids_a, ids_a),
                              sample_metadata = x$sample_metadata),
    "batch adjustment")
})

test_that("stability refits carry adjust_batches", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  skip_if_not_installed("sva")
  x <- .toy()
  fit <- canonicalize_programs(fit_celltype_programs(x, adjust_batches = "batch",
                                                     max_factors = 2, seed = 1))
  st <- assess_program_stability(fit, n_boot = 2, resample_frac = 0.8,
                                 seed = 1, progress = FALSE)
  expect_gt(nrow(st$stability$per_program), 0)
})

test_that("span_alignment recovers identical, subset, and orthogonal spans", {
  set.seed(5)
  W <- matrix(rnorm(200 * 3), 200, 3)
  a <- span_alignment(W, W)
  expect_equal(a$principal_cosines, rep(1, 3), tolerance = 1e-8)
  expect_equal(a$mean_overlap, 1)
  # subset span: the 2 leading cosines are 1 (fully contained), third is partial
  b <- span_alignment(W, W[, 1:2, drop = FALSE])
  expect_equal(b$principal_cosines[1:2], rep(1, 2), tolerance = 1e-8)
  expect_equal(b$dim, c(3, 2))
  # orthogonal spans
  Wx <- W; Wx[, 1] <- W[, 2]; Wx[, 2] <- W[, 1]  # same span, permuted: still 1
  expect_equal(span_alignment(W, Wx)$mean_overlap, 1, tolerance = 1e-8)
  Wz <- matrix(rnorm(200 * 2), 200, 2)  # random span: overlap well below 1
  expect_lt(span_alignment(W, Wz)$mean_overlap, 0.9)
  # row-name matching: permuted rows realigned by name
  rownames(W) <- paste0("g", seq_len(200))
  Wp <- W[sample(200), ]
  expect_equal(span_alignment(W, Wp)$mean_overlap, 1, tolerance = 1e-8)
  expect_error(span_alignment(W, W[1:50, ]), "same number")
})
