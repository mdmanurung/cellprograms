test_that("CPUB-004 historical node counts separate events and failure accounting", {
  # The portable oracle must be installed: absence is an error, never a skip.
  path <- system.file("validation", "analyze_corrected.R",
                      package = "cellprograms", mustWork = TRUE)
  benchmark <- testthat::test_path("..", "..", "benchmarks", "cmp_calibration",
                                  "analyze_corrected.R")
  if (file.exists(benchmark)) {
    expect_identical(readLines(path), readLines(benchmark))
  }
  env <- new.env(parent = globalenv()); sys.source(path, env)
  requested <- data.frame(scenario = "mixed", seed = 1:4)
  r <- data.frame(scenario = "mixed", seed = 1:3, status = c("ok", "ok", "error"),
    null_n = c(100, 100, NA), null_n_sig = c(0, 1, NA),
    plant_n = c(2, 2, NA), plant_n_sig = c(2, 2, NA))
  out <- env$summarize_node_counts(r, requested)
  expect_equal(out$n_requested, 4)
  expect_equal(out$n_callable, 2)
  expect_equal(out$M06_conditional_node_family, .5)
  expect_equal(out$M07_pooled_null_node_FPR, .005)
  expect_equal(out$M06_failure_inclusive_upper_rate, .75)
  expect_equal(out$M06_conditional_upper95, sqrt(.95), tolerance = 1e-10)
  expect_true(is.na(out$edge_or_set_FWER))
  expect_false(out$calibration_certified)
  expect_error(env$summarize_node_counts(r, requested[c(1, 1), ]), "unique")
})
