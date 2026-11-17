## load_package.R -- make the package (internal helpers included) available to a benchmark script.
## The scripts used to source six files of R/ by name. The package has since grown helpers in other
## files (covariate-replay.R, preprocess-freeze.R, ...), so a fixed list no longer works.
load_package <- function(root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE, export_all = TRUE)
  } else {
    for (f in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
  }
  invisible(TRUE)
}
