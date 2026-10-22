# CI policy oracles. These run from the repository, not the installed package.
policy <- normalizePath(".github/check-core-notes.R")
tmp <- tempfile("note-policy-")
dir.create(tmp)
old <- setwd(tmp)
tryCatch({
  dir.create("check")
  check <- function(text, should_fail) {
    writeLines(text, "check/00check.log")
    failed <- tryCatch({ source(policy, local = new.env()); FALSE },
                       error = function(e) TRUE)
    stopifnot(identical(failed, should_fail))
  }
  check(c("* checking package dependencies ... NOTE",
          "Packages suggested but not available for checking: 'Optional'",
          "* checking examples ... OK", "Status: 1 NOTE"), FALSE)
  check(c("* checking package dependencies ... OK",
          "* checking examples ... OK", "Status: OK"), FALSE)
  check(c("* checking R code ... NOTE", "Unregistered native call",
          "* checking examples ... OK", "Status: 1 NOTE"), TRUE)
  check(c("* checking examples ... WARNING", "A warning", "Status: 1 WARNING"), TRUE)
  cat("Four core-note policy oracles passed.\n")
}, finally = {
  setwd(old)
  unlink(tmp, recursive = TRUE)
})
