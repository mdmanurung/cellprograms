# Regenerates fixtures/sharing_perf_golden.rds for test-sharing-perf-regression.R.
# Run from the package root: Rscript tests/testthat/fixtures/make-sharing-perf-golden.R
# The committed file was captured from release/0.3.0 @ 0914144 (loop with cor() on
# permuted copies; release/0.3.0 and the current loop agree). Do not regenerate it to make a failing test pass.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("tests/testthat/helper-sharing-perf.R")
golden <- list(
  continuous = .match_programs_loadings(.sharing_perf_fit(), n_perm = 199L, seed = 1L),
  ties = .match_programs_loadings(.sharing_tie_fit(), n_perm = 99L, seed = 1L)
)
saveRDS(golden, "tests/testthat/fixtures/sharing_perf_golden.rds", version = 2)
