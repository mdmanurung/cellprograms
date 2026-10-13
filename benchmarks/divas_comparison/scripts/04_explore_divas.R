# 04_explore_divas.R -- inspect the structure of the precomputed DIVAS results
# (run once to design 04_extract_divas.R; not part of the final pipeline)
#
# Usage: Rscript 04_explore_divas.R <data_dir>

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"

obj <- readRDS(file.path(data_dir, "divas_results_combined_6block_renamed.rds"))
cat("class:", class(obj), "\n")
cat("top-level names:\n")
print(names(obj))

explore <- function(x, name, depth = 1) {
  if (depth > 3) return(invisible())
  if (is.list(x)) {
    cat(strrep(" ", (depth - 1) * 2), name, " : list(", length(x), "): ",
        paste(utils::head(names(x), 12), collapse = ", "),
        if (length(x) > 12) " ..." else "", "\n", sep = "")
    for (nm in utils::head(names(x), 8)) {
      explore(x[[nm]], nm, depth + 1)
    }
  } else if (is.matrix(x) || is.data.frame(x)) {
    cat(strrep(" ", (depth - 1) * 2), name, " : ", class(x)[1], " dim=",
        paste(dim(x), collapse = "x"),
        " dimnames[1][1:4]=", paste(utils::head(rownames(x), 4), collapse = ","),
        " dimnames[2][1:6]=", paste(utils::head(colnames(x), 6), collapse = ","),
        "\n", sep = "")
  } else if (is.vector(x) && length(x) <= 8) {
    cat(strrep(" ", (depth - 1) * 2), name, " : ", paste(x, collapse = ","), "\n", sep = "")
  } else {
    cat(strrep(" ", (depth - 1) * 2), name, " : ", class(x)[1],
        " len=", length(x), "\n", sep = "")
  }
}

for (nm in names(obj)) explore(obj[[nm]], nm)
