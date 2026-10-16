## Fit donor-level programs on pseudobulks via cellprograms (R).
##
## Reads pseudobulk.h5 (rhdf5), log1p-CP10K normalizes each donor-ct row, fits
## each cell type SEPARATELY (per-ct tryCatch so one degenerate cell type cannot
## kill the run; failures logged), merges successful fits into one
## cell_program_fit object, saves fit RDS + summary CSV + failures CSV.
##
## Usage: Rscript fit_pilot.R <pseudobulk_dir> <cohort> <out_rds> [cts_file]
suppressPackageStartupMessages({
  library(cellprograms)
  library(rhdf5)
})

args <- commandArgs(TRUE)
pb_dir <- args[1]
cohort <- args[2]
out_rds <- args[3]
cts_filter <- if (length(args) >= 4 && nzchar(args[4])) {
  readLines(args[4])
} else NULL

.safe <- function(s) {  # mirror Python safe_name
  z <- gsub("[^A-Za-z0-9_.-]+", "_", s)
  z <- sub("^_+", "", z)
  z <- sub("_+$", "", z)
  substr(z, 1, 80)
}

h5file <- file.path(pb_dir, "pseudobulk.h5")
genes <- h5read(h5file, "genes")
cts_all <- h5read(h5file, "cell_types")
cts <- if (is.null(cts_filter)) cts_all else intersect(cts_all, cts_filter)
stopifnot(length(cts) > 0)

meta <- read.csv(file.path(pb_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
meta$observation_id <- meta$donor_id

pseudobulk <- list()
for (ct in cts) {
  g <- h5read(h5file, paste0("/", .safe(ct), "/counts"))
  donors_ct <- h5read(h5file, paste0("/", .safe(ct), "/donors"))
  if (nrow(g) == length(genes) && ncol(g) != length(genes)) g <- t(g)  # rhdf5 C-order
  rownames(g) <- donors_ct
  colnames(g) <- genes
  lib <- rowSums(g)
  lib[lib == 0] <- 1
  Y <- log1p(g / lib * 1e4)
  # reduce to the top-2000 variable genes NOW (matches .select_features
  # "variable" ranking on the same log-normalized values) so big cohorts
  # (981 donors x 36k genes x 58 cts) fit in memory; the fit re-selects all.
  v <- apply(Y, 2, stats::var)
  k <- min(2000L, length(v))
  pseudobulk[[ct]] <- Y[, names(sort(v, decreasing = TRUE))[seq_len(k)], drop = FALSE]
  rm(g, Y, v); invisible(gc())
}
rhdf5::H5close()

x_full <- as_cell_program_data(pseudobulk, meta)

fit_one <- function(ct) {
  x_ct <- x_full
  x_ct$pseudobulk <- pseudobulk[ct]
  x_ct$cell_types <- ct
  x_ct$observation_ids <- rownames(pseudobulk[[ct]])
  f <- fit_celltype_programs(
    x_ct,
    loading_prior = "point_laplace",
    max_factors = 10,
    features = "variable",
    n_variable_genes = 2000,
    var_type = 1,
    feature_scale = "logcounts",
    seed = 20261016
  )
  canonicalize_programs(f)
}

fits <- list(); failures <- list()
for (ct in cts) {
  f <- tryCatch(fit_one(ct), error = function(e) {
    message("FIT FAILED [", ct, "]: ", conditionMessage(e))
    NULL
  })
  if (is.null(f)) {
    failures[[length(failures) + 1]] <- data.frame(cell_type = ct)
  } else {
    fits[[ct]] <- f
  }
}
if (!length(fits)) stop("all cell-type fits failed")
cat(sprintf("fitted %d/%d cell types (%d failed)\n", length(fits), length(cts), length(failures)))

# merge per-ct fits into one object (template = first successful fit)
merged <- fits[[1]]
if (length(fits) > 1) {
  for (ct in names(fits)[-1]) {
    f <- fits[[ct]]
    for (nm in c("fits", "programs", "scores", "loadings", "canonicalization")) {
      if (!is.null(f[[nm]])) merged[[nm]][[ct]] <- f[[nm]][[ct]]
    }
    merged$cell_types <- c(merged$cell_types, ct)
  }
}
# global fields from the full data object
merged$pseudobulk <- x_full$pseudobulk
merged$sample_metadata <- x_full$sample_metadata
merged$observation_ids <- x_full$observation_ids
merged$shared_genes <- x_full$shared_genes
merged$dims <- x_full$dims

saveRDS(merged, out_rds)

summ <- program_summary(merged)
write.csv(summ, sub("\\.rds$", "_summary.csv", out_rds), row.names = FALSE)
if (length(failures)) {
  write.csv(do.call(rbind, failures), sub("\\.rds$", "_failures.csv", out_rds), row.names = FALSE)
}
cat("fit saved:", out_rds, "\n")
print(summ)
