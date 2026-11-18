# Fixture for the .match_programs_loadings() regression guard
# (test-sharing-perf-regression.R). Loadings are planted; no flashier needed.
# A_1 = B_1 (shared), C_2 and D_3 weakly related, everything else unrelated, so the
# permutation p-values take a spread of values rather than all sitting at 1/(n_perm + 1).
.sharing_perf_fit <- function(g = 400, k = 4, seed = 4, cts = c("A", "B", "C", "D", "E")) {
  set.seed(seed)
  genes <- sprintf("G%03d", seq_len(g))
  shared <- rnorm(g)
  weak <- rnorm(g)
  loadings <- fits <- canon <- list()
  for (ct in cts) {
    W <- matrix(rnorm(g * k, sd = 1), g, k, dimnames = list(genes, paste0(ct, "_", seq_len(k))))
    if (ct %in% c("A", "B")) W[, 1] <- shared + rnorm(g, sd = 0.2)
    if (ct == "C") W[, 2] <- W[, 2] + 0.5 * weak
    if (ct == "D") W[, 3] <- W[, 3] + 0.5 * weak
    loadings[[ct]] <- W
    fits[[ct]] <- list(flash = list(F_pm = W, F_psd = matrix(0.3, g, k, dimnames = dimnames(W))),
                       genes = genes)
    canon[[ct]] <- list(scale_multiplier = rep(1, k), sign_multiplier = rep(1, k))
  }
  structure(list(loadings = loadings, fits = fits, canonicalization = canon, cell_types = cts),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

# Discrete loadings make exact ties |null| == |obs| common. Different arithmetic for the
# same correlation (cor() vs crossprod of scaled columns) breaks such ties differently,
# so this fixture guards the p-values against a change in how the null is computed.
.sharing_tie_fit <- function(g = 30, k = 2, seed = 5, cts = c("A", "B")) {
  set.seed(seed)
  genes <- sprintf("G%02d", seq_len(g))
  loadings <- fits <- canon <- list()
  for (ct in cts) {
    W <- matrix(sample(c(-1, 0, 1, 2), g * k, TRUE), g, k,
                dimnames = list(genes, paste0(ct, "_", seq_len(k))))
    loadings[[ct]] <- W
    fits[[ct]] <- list(flash = list(F_pm = W, F_psd = matrix(0.3, g, k, dimnames = dimnames(W))),
                       genes = genes)
    canon[[ct]] <- list(scale_multiplier = rep(1, k), sign_multiplier = rep(1, k))
  }
  structure(list(loadings = loadings, fits = fits, canonicalization = canon, cell_types = cts),
            class = c("cell_program_fit", "cell_program_data", "list"))
}
