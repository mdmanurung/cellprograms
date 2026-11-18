# Purpose: time the permutation loops that dominate find_coordinated_programs() and
#   program_sharing(), so a change can be compared against an older commit.
# Run: Rscript benchmarks/performance/bench_cmp.R [PACKAGE_ROOT]
#   PACKAGE_ROOT defaults to "."; point it at a worktree of an older commit to time
#   that code. The fixtures always come from the current tree's tests/testthat/.
# Inputs:
#   tests/testthat/helper-toy.R           .cmp_visit_toy()
#   tests/testthat/helper-sharing-perf.R  .sharing_perf_fit()
# Outputs:
#   one line per case on stdout (best of N_REP wall-clock seconds)

# Parameters -----------------------------------------------------------------------------
N_REP <- 3          # repeats per case; the minimum is reported
N_DONORS <- 200     # donors in the coordination cases
K_PROGRAMS <- 5     # programs per cell type
N_CELL_TYPES <- c(5, 10, 20)   # coordination scale
N_PERM_CMP <- 99    # permutations (the default, 199, scales linearly)

# Setup ----------------------------------------------------------------------------------
root <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(root)) root <- "."
helpers <- normalizePath("tests/testthat")
suppressMessages(pkgload::load_all(root, quiet = TRUE))
source(file.path(helpers, "helper-toy.R"))
source(file.path(helpers, "helper-sharing-perf.R"))
# re-evaluates the expression each repeat (a bare promise would be forced only once)
best <- function(expr) {
  call <- substitute(expr)
  env <- parent.frame()
  min(replicate(N_REP, system.time(eval(call, env))[["elapsed"]]))
}

# find_coordinated_programs --------------------------------------------------------------
for (n_ct in N_CELL_TYPES) {
  fit <- .cmp_visit_toy(n_pat = N_DONORS, visits = "T1", k = K_PROGRAMS,
                        cts = paste0("ct", seq_len(n_ct)))
  fit$sample_metadata$site <- rep(c("a", "b"), length.out = nrow(fit$sample_metadata))
  secs <- best(find_coordinated_programs(fit, strata = "site", n_perm = N_PERM_CMP,
                                         ambient = "none"))
  cat(sprintf("find_coordinated_programs  %2d cell types, %d donors, n_perm = %d: %6.2f s\n",
              n_ct, N_DONORS, N_PERM_CMP, secs))
}

# program_sharing's loading match --------------------------------------------------------
fit <- .sharing_perf_fit(g = 2000, k = K_PROGRAMS, cts = paste0("ct", 1:10))
secs <- best(.match_programs_loadings(fit, n_perm = 199L, seed = 1L))
cat(sprintf(".match_programs_loadings   10 cell types, 2000 genes, n_perm = 199: %6.2f s\n", secs))
