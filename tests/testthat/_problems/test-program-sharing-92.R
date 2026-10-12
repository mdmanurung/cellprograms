# Extracted from test-program-sharing.R:92

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
  common <- rnorm(g)   # shared across cell types: plants CROSS-CELL-TYPE null correlation
  for (ct in cts) {
    K <- length(eff[[ct]])
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
