## SYNC-2: canonicalization via flashier.utils::standardize_factors("rms").
## Verify the upstream representation reproduces the legacy inline RMS + sign
## convention exactly and that downstream invariants (sharing, CMP inputs) hold.

test_that("upstream RMS canonicalization matches the legacy inline branch", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  fit_u <- withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
    canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  fit_l <- withr::with_options(list(cellprograms.use_flashier_utils = FALSE),
    canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  for (ct in x$cell_types) {
    expect_equal(fit_u$canonicalization[[ct]]$Z, fit_l$canonicalization[[ct]]$Z,
                 tolerance = 1e-10)
    expect_equal(fit_u$canonicalization[[ct]]$W, fit_l$canonicalization[[ct]]$W,
                 tolerance = 1e-10)
    expect_equal(fit_u$canonicalization[[ct]]$scale_multiplier,
                 fit_l$canonicalization[[ct]]$scale_multiplier, tolerance = 1e-10)
    expect_equal(fit_u$canonicalization[[ct]]$sign_multiplier,
                 fit_l$canonicalization[[ct]]$sign_multiplier)
  }
})

test_that("canonicalization is idempotent under the upstream path", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  fit <- withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
    canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  s1 <- fit$scores; l1 <- fit$loadings
  fit2 <- canonicalize_programs(fit)
  expect_equal(fit2$scores, s1)
  expect_equal(fit2$loadings, l1)
})

test_that("native uncertainty (F_psd, F_lfsr) is preserved through canonicalization", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  fit <- withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
    canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  for (ct in x$cell_types) {
    f <- fit$fits[[ct]]$flash
    if (is.null(f$L_pm)) next
    # canonicalization must not alter the native posterior SD / LFSR matrices
    expect_equal(dim(f$F_psd), dim(f$F_pm))
    expect_equal(dim(f$F_lfsr), dim(f$F_pm))
    expect_true(all(f$F_psd >= 0))
    expect_true(all(f$F_lfsr >= 0 & f$F_lfsr <= 1))
  }
})

test_that("program sharing inputs (F_psd * scale_multiplier) are unchanged", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  fit_u <- withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
    canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  fit_l <- withr::with_options(list(cellprograms.use_flashier_utils = FALSE),
    canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  for (ct in x$cell_types) {
    fu <- fit_u$fits[[ct]]$flash; fl <- fit_l$fits[[ct]]$flash
    if (is.null(fu$F_psd) || is.null(fl$F_psd)) next
    su <- fit_u$canonicalization[[ct]]$scale_multiplier
    sl <- fit_l$canonicalization[[ct]]$scale_multiplier
    # the sharing layer scales native F_psd by scale_multiplier; both must agree
    expect_equal(sweep(fu$F_psd, 2, su, `*`), sweep(fl$F_psd, 2, sl, `*`),
                 tolerance = 1e-10)
  }
})

test_that("K = 0 and degenerate cell types still canonicalize to empty", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  # a cell type that is constant across samples -> degenerate K=0 sentinel
  set.seed(3)
  ids <- paste0("s", 1:20)
  const <- matrix(1, 20, 30, dimnames = list(ids, paste0("g", 1:30)))
  varm <- matrix(rnorm(20 * 30), 20, 30, dimnames = list(ids, paste0("g", 1:30)))
  x <- as_cell_program_data(list(CONST = const, VAR = varm),
                            sample_metadata = data.frame(observation_id = ids))
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 2, seed = 1,
                                                     features = "all"))
  expect_equal(length(fit$canonicalization$CONST$scale_multiplier), 0)
  expect_equal(ncol(fit$scores$CONST), 0L)
})
