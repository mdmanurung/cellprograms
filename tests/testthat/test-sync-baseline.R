## SYNC-0 baseline: freeze the ordinary-fit output contract before the
## flashier.utils migration. These tests record current behavior; SYNC-1B/2 must
## keep them passing (numerical + scientific non-regression).

test_that("ordinary fit exposes the documented output contract", {
  fit <- sync_baseline_fit()
  expect_s3_class(fit, "cell_program_fit")
  for (field in c("pseudobulk", "sample_metadata", "cell_types", "fits",
                  "preprocessing", "canonicalization", "scores", "loadings")) {
    expect_true(field %in% names(fit), info = paste("missing field:", field))
  }
  for (ct in fit$cell_types) {
    expect_true(all(c("flash", "Y_used", "genes", "fallback") %in% names(fit$fits[[ct]])))
    expect_true(all(c("Z", "W", "scale_multiplier", "sign_multiplier") %in%
                      names(fit$canonicalization[[ct]])))
  }
})

test_that("canonicalization preserves reconstruction exactly", {
  fit <- sync_baseline_fit()
  for (ct in fit$cell_types) {
    f <- fit$fits[[ct]]$flash
    if (is.null(f$L_pm) || is.null(f$F_pm)) next
    recon_native <- tcrossprod(as.matrix(f$L_pm), as.matrix(f$F_pm))
    recon_canon <- tcrossprod(fit$canonicalization[[ct]]$Z, fit$canonicalization[[ct]]$W)
    expect_equal(recon_canon, recon_native, tolerance = 1e-8)
  }
})

test_that("RMS canonicalization uses unit-RMS scores and compensating loadings", {
  fit <- sync_baseline_fit()
  for (ct in fit$cell_types) {
    Z <- fit$canonicalization[[ct]]$Z
    if (is.null(Z) || ncol(Z) == 0L) next
    rms <- sqrt(colMeans(Z^2))
    expect_equal(unname(rms), rep(1, ncol(Z)), tolerance = 1e-8)
  }
})

test_that("sign convention anchors on the max-|loading| gene, deterministic", {
  fit <- sync_baseline_fit()
  for (ct in fit$cell_types) {
    W <- fit$canonicalization[[ct]]$W
    if (is.null(W) || ncol(W) == 0L) next
    for (k in seq_len(ncol(W))) {
      anchor <- which.max(abs(W[, k]))
      expect_gt(W[anchor, k], 0)
    }
  }
})

test_that("program IDs follow the <cell_type>_<k> convention", {
  fit <- sync_baseline_fit()
  for (ct in fit$cell_types) {
    K <- ncol(fit$scores[[ct]])
    expect_equal(colnames(fit$scores[[ct]]), paste0(ct, "_", seq_len(K)))
    expect_equal(colnames(fit$loadings[[ct]]), paste0(ct, "_", seq_len(K)))
  }
})

test_that("baseline fingerprint is stable across repeated fits (seeded)", {
  f1 <- sync_fingerprint(sync_baseline_fit())
  f2 <- sync_fingerprint(sync_baseline_fit())
  expect_equal(f1$recon_max_abs_err, 0, tolerance = 1e-8)
  expect_equal(f1$scale_multiplier, f2$scale_multiplier, tolerance = 1e-10)
  expect_equal(f1$sign_multiplier, f2$sign_multiplier)
  expect_equal(f1$n_programs, f2$n_programs)
  expect_equal(f1$fallback, f2$fallback)
})
