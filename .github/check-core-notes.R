# Development ordinary-core policy, not an --as-cran release certificate.
# Only an unavailable optional-Suggests dependency note is allowed.
logs <- list.files("check", pattern = "^00check\\.log$", recursive = TRUE,
                   full.names = TRUE)
if (!length(logs)) {
  logs <- list.files(".", pattern = "^00check\\.log$", recursive = TRUE,
                     full.names = TRUE)
}
if (length(logs) != 1L) stop("Expected exactly one retained R CMD check log.")
lines <- readLines(logs[[1]], warn = FALSE)
starts <- grep("^\\* checking .* \\.\\.\\. ", lines)
for (i in seq_along(starts)) {
  from <- starts[[i]]
  to <- if (i < length(starts)) starts[[i + 1L]] - 1L else length(lines)
  block <- paste(lines[from:to], collapse = "\n")
  if (grepl("\\.\\.\\. (WARNING|ERROR)$", lines[[from]])) stop(block)
  if (grepl("\\.\\.\\. NOTE$", lines[[from]])) {
    allowed <- grepl("^\\* checking package dependencies", lines[[from]]) &&
      grepl("suggested.*not available|Suggests.*not available", block,
            ignore.case = TRUE)
    if (!allowed) stop("Unexpected core-check note:\n", block)
  }
}
cat("Development core notes: only explicitly unavailable optional Suggests allowed.\n")
