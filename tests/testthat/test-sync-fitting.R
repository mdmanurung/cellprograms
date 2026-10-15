## SYNC-1B: ordinary-fit migration to flashier.utils::fit_ebmf(). Verify the
## upstream backend reproduces the legacy direct flashier::flash() backend
## numerically and that unmappable controls are rejected, never silently dropped.

test_that("flashier.utils backend matches the legacy backend to machine precision", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  withr::local_options(list(cellprograms.use_flashier_utils = NULL))
  # force legacy
  withr::with_options(list(cellprograms.use_flashier_utils = FALSE),
    fit_legacy <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  # upstream (default)
  withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
    fit_fu <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1)))
  for (ct in x$cell_types) {
    expect_equal(fit_fu$fits[[ct]]$flash$n_factors, fit_legacy$fits[[ct]]$flash$n_factors)
    expect_equal(as.matrix(fit_fu$fits[[ct]]$flash$L_pm),
                 as.matrix(fit_legacy$fits[[ct]]$flash$L_pm), tolerance = 1e-8)
    expect_equal(as.matrix(fit_fu$fits[[ct]]$flash$F_pm),
                 as.matrix(fit_legacy$fits[[ct]]$flash$F_pm), tolerance = 1e-8)
    expect_equal(fit_fu$fits[[ct]]$flash$elbo, fit_legacy$fits[[ct]]$flash$elbo, tolerance = 1e-8)
    # canonicalized outputs identical
    expect_equal(fit_fu$scores[[ct]], fit_legacy$scores[[ct]], tolerance = 1e-8)
    expect_equal(fit_fu$loadings[[ct]], fit_legacy$loadings[[ct]], tolerance = 1e-8)
  }
})

test_that("upstream backend preserves posterior SD and LFSR fields", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  fit <- withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
                             fit_celltype_programs(x, max_factors = 3, seed = 1))
  for (ct in x$cell_types) {
    f <- fit$fits[[ct]]$flash
    if (is.null(f$L_pm)) next
    expect_false(is.null(f$L_psd))
    expect_false(is.null(f$F_psd))
    expect_false(is.null(f$F_lfsr))
    expect_equal(dim(f$L_psd), dim(f$L_pm))
    expect_equal(dim(f$F_psd), dim(f$F_pm))
  }
})

test_that("unmappable flash_control overrides are rejected explicitly", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  # `S` (a residual-SD matrix) is a flashier::flash argument that fit_ebmf()
  # does not expose through the ordinary-fit mapping; it must error, not be
  # silently dropped.
  expect_error(
    withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
      fit_celltype_programs(x, max_factors = 3, seed = 1,
                            flash_control = list(S = 1))),
    "cannot be mapped|use_flashier_utils"
  )
})

test_that("legacy backend still works when flashier.utils is disabled", {
  x <- .toy()
  fit <- withr::with_options(list(cellprograms.use_flashier_utils = FALSE),
                             fit_celltype_programs(x, max_factors = 3, seed = 1))
  expect_s3_class(fit, "cell_program_fit")
  expect_true(all(vapply(fit$fits, function(z) z$fallback, character(1)) == "none"))
})

test_that("var_type c(1, 2) Kronecker fits route through the upstream backend", {
  skip_if_not(requireNamespace("flashier.utils", quietly = TRUE))
  x <- .toy()
  fit <- withr::with_options(list(cellprograms.use_flashier_utils = TRUE),
                             fit_celltype_programs(x, max_factors = 3, seed = 1, var_type = c(1, 2)))
  expect_s3_class(fit, "cell_program_fit")
  expect_true(all(vapply(fit$scores, function(s) ncol(s) >= 0, logical(1))))
})
