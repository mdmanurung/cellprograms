## Shared fixtures, memoized to keep the suite fast.
.test_cache <- new.env()

get_sim_fit <- function(scenario = 4L, N = 100L, G = 300L, seed = 42L) {
  key <- paste(scenario, N, G, seed, sep = "-")
  if (!exists(key, envir = .test_cache)) {
    sim <- simulate_programs(scenario = scenario, N = N, G = G, signal_sd = 3,
                             sparsity = 0.05, noise_sd = 1, seed = seed)
    fit <- canonicalize_programs(fit_celltype_programs(sim$x, max_factors = 10L))
    assign(key, list(sim = sim, fit = fit), envir = .test_cache)
  }
  get(key, envir = .test_cache)
}

## max |cor| between a cell type's programs and a true activity vector
## (simulate_programs observation ids are "obs<i>", positional in z)
best_cor <- function(fit, ct, z) {
  S <- fit$scores[[ct]]
  if (is.null(S) || ncol(S) == 0L) return(NA_real_)
  pos <- as.integer(sub("^obs", "", rownames(S)))
  max(abs(stats::cor(S, z[pos])))
}
