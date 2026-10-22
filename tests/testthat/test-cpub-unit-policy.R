test_that("unknown reference units require refit before any native projection", {
  # The guard is before native inference. No unit labels are inferred or changed.
  x <- as_cell_program_data(list(A = matrix(
    seq_len(48), 8, 6, dimnames = list(paste0("d", 1:8), paste0("g", 1:6)))))
  fit <- x
  class(fit) <- c("cell_program_fit", class(x))
  fit$scores <- list(A = matrix(1, 8, 1))
  fit$canonicalization <- list(A = list(scale_multiplier = 1, sign_multiplier = 1))
  ct <- setNames(rep("A", 8), rownames(x$pseudobulk$A))
  for (units in list(NULL, "unknown", NA_character_, "", character())) {
    fit$preprocessing <- list(feature_scale = units)
    original <- serialize(fit, NULL)
    for (query in list(NULL, "synthetic_linear")) {
      expect_error(project_programs_to_cells(fit, x$pseudobulk$A, ct,
                                             feature_scale = query),
                   "refit with truthful")
    }
    expect_identical(serialize(fit, NULL), original)
  }
})

test_that("warning inventory assertions fail on unexpected or missing warnings", {
  # Test the test helper itself so warning counts cannot be reduced silently.
  ok <- testthat::capture_expectation(
    .expect_warning_messages(warning("expected"), "expected"))
  expect_s3_class(ok, "expectation_success")
  extra <- testthat::capture_expectation(
    .expect_warning_messages({ warning("expected"); warning("extra") }, "expected"))
  expect_s3_class(extra, "expectation_failure")
  absent <- testthat::capture_expectation(.expect_warning_messages(1, "expected"))
  expect_s3_class(absent, "expectation_failure")
})
