# Every captured warning must match the explicitly expected complete message set.
# Additional, missing or changed warnings fail; this is not blanket suppression.
.expect_warning_messages <- function(code, expected) {
  seen <- character()
  value <- withCallingHandlers(force(code), warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  testthat::expect_identical(sort(seen), sort(expected))
  value
}
