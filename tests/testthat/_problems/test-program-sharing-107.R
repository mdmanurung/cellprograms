# Extracted from test-program-sharing.R:107

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "cellprograms", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.sharing_toy <- function(g = 300, seed = 1, null_rho = 0) {
  set.seed(seed)
  genes <- sprintf("G%03d", seq_len(g))
  cts <- c("A", "B", "C")
  eff <- list(
    A = list(c(rnorm(40, sd = 2), rep(0, g - 40)),                 # A_1: shared
             c(rep(0, 60), rnorm(40, sd = 2), rep(0, g - 100))),   # A_2: private
    B = list(NULL),
    C = list(c(rep(0, 120), rnorm(40, sd = 2), rep(0, g - 160)))   # C_1: private
  )
  eff$B[[1]] <- eff$A[[1]]                                         # identical shared program
  loadings <- fits <- canon <- list()
  for (ct in cts) {
    K <- length(eff[[ct]])
    common <- rnorm(g)
    W <- matrix(rnorm(g * K, sd = 0.05), g, K,
                dimnames = list(genes, paste0(ct, "_", seq_len(K))))
    for (k in seq_len(K)) {
      W[, k] <- W[, k] + eff[[ct]][[k]] + sqrt(null_rho) * 0.05 * common
    }
    loadings[[ct]] <- W
    fits[[ct]] <- list(
      flash = list(F_pm = W, F_psd = matrix(0.3, g, K, dimnames = dimnames(W))),
      genes = genes
    )
    canon[[ct]] <- list(scale_multiplier = rep(1, K), sign_multiplier = rep(1, K))
  }
  structure(list(loadings = loadings, fits = fits, canonicalization = canon,
                 cell_types = cts),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

# test -------------------------------------------------------------------------
skip_if_not_installed("mashr")
fit <- .sharing_toy()
res <- program_sharing(fit, n_perm = 199, seed = 1)
expect_s3_class(res, "program_sharing")
s <- res$summary
expect_equal(s$classification[s$program_id == "A_1"], "shared")
expect_equal(s$classification[s$program_id == "B_1"], "shared")
expect_equal(s$classification[s$program_id == "A_2"], "private")
expect_equal(s$classification[s$program_id == "C_1"], "private")
expect_equal(s$evidence[s$program_id == "C_1"], "no_significant_match")
r <- res$results$group1
expect_equal(dim(r$posterior_mean), c(300L, 2L))
expect_equal(dim(r$lfsr), c(300L, 2L))
expect_gt(r$p_all_shared, 0.5)
expect_gt(sum(r$pattern_weights), 0.99)
expect_true("equal_effects" %in% names(r$pattern_weights))
forced <- list(data.frame(cell_type = c("A", "C"), program_id = c("A_2", "C_1")))
res2 <- program_sharing(fit, matches = forced)
expect_equal(res2$results$group1$classification, "private")
expect_equal(res2$results$group1$frac_top_private, 1)
