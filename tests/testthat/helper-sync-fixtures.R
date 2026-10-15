## SYNC-0 baseline fixtures: freeze the ordinary-fit output contract so the
## flashier.utils migration (SYNC-1B/2) can be checked for numerical and
## scientific non-regression against a recorded reference.

## A small deterministic ordinary fit (no covariates) used as the parity
## reference. max_factors = 3 keeps it fast; seed fixed for reproducibility.
sync_baseline_fit <- function() {
  x <- .toy()
  fit <- fit_celltype_programs(x, max_factors = 3, seed = 1)
  canonicalize_programs(fit)
}

## Scalar + structural fingerprint of a fit's ordinary-fit contract: the fields
## the migration must preserve. Returns a flat list suitable for expect_equal.
sync_fingerprint <- function(fit) {
  ct <- fit$cell_types[1]
  f <- fit$fits[[ct]]$flash
  list(
    cell_types = fit$cell_types,
    n_programs = vapply(fit$scores, ncol, 1L),
    score_dim = lapply(fit$scores, dim),
    loading_dim = lapply(fit$loadings, dim),
    fallback = vapply(fit$fits, function(z) z$fallback, character(1)),
    scale_multiplier = fit$canonicalization[[ct]]$scale_multiplier,
    sign_multiplier = fit$canonicalization[[ct]]$sign_multiplier,
    # reconstruction identity: Z %*% t(W) must equal the native L %*% t(F)
    recon_max_abs_err = max(abs(
      tcrossprod(fit$canonicalization[[ct]]$Z, fit$canonicalization[[ct]]$W) -
        tcrossprod(as.matrix(f$L_pm), as.matrix(f$F_pm))
    )),
    scores_colnames = lapply(fit$scores, colnames)
  )
}
