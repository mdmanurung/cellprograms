test_that("maintainer and installed public support metadata are usable", {
  d <- utils::packageDescription("cellprograms")
  expect_false(grepl("example\\.(org|com)", d$`Authors@R`))
  expect_match(d$License, "MIT")
  support <- system.file("validation", "CORE_SUPPORT.md", package = "cellprograms",
                         mustWork = TRUE)
  text <- paste(readLines(support), collapse = "\n")
  expect_match(text, "Unknown-unit reference projection")
  expect_match(text, "Optional")
})
