test_that("CMP printing states nominal adjustment rather than certified FWER", {
  res <- structure(list(
    cmps = .cmp_empty(),
    settings = list(alpha = .05, mode = "marginal", ambient = "none")
  ), class = "cmp_result")
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "nominal node threshold 0.05", fixed = TRUE)
  expect_match(out, "Fixed-input permutation-max adjustment", fixed = TRUE)
  expect_match(out, "CMP-set error control not certified", fixed = TRUE)
  expect_false(grepl("at FWER", out, fixed = TRUE))
  invisible(capture.output(returned <- withVisible(print(res))))
  expect_false(returned$visible)
  expect_identical(returned$value, res)
})

test_that("historical cap notification does not prescribe a general cap", {
  old <- options(cellprograms.cmp_cap_message_shown = NULL)
  on.exit(options(old))
  fit <- list(preprocessing = list(max_factors = 2),
              scores = list(A = matrix(0, 4, 2)))
  msg <- capture.output(
    withCallingHandlers(.cmp_cap_message(fit),
      message = function(m) {
        cat(conditionMessage(m)); invokeRestart("muffleMessage")
      }))
  text <- paste(msg, collapse = "\n")
  expect_match(text, "summary-only evidence", fixed = TRUE)
  expect_match(text, "not a general cap recommendation", fixed = TRUE)
  expect_true(isTRUE(getOption("cellprograms.cmp_cap_message_shown")))
})
