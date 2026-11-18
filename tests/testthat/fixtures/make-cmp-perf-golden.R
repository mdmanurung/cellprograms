# Regenerates fixtures/cmp_perf_golden.rds for test-cmp-perf-regression.R.
# Run from the package root on a tree whose find_coordinated_programs() you trust:
#   Rscript tests/testthat/fixtures/make-cmp-perf-golden.R
# The committed file was captured from release/0.3.0 @ 0914144, before the
# permutation loop was rewritten. Do not regenerate it to make a failing test pass.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("tests/testthat/helper-cmp-perf.R")
golden <- lapply(.perf_cases(), \(case) .perf_pin(.perf_run(case)))
golden$stability <- .perf_stability_pin()
saveRDS(golden, "tests/testthat/fixtures/cmp_perf_golden.rds", version = 2)
