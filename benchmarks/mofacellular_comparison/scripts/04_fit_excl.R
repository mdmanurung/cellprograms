## 04_fit_excl.R -- sensitivity refit: exclude_artifacts = TRUE (all 4 sets).
## Identical recipe to 01_fit_masked.R but with the full artifact preset
## (vdj + genotype_marker + sex + canonical_marker) instead of the 274-gene
## vdj + genotype_marker list. Exclusion still happens before HVG selection.
##
## Usage: Rscript 04_fit_excl.R [fit_rds] [out_rds]
##   defaults: results/stephenson_K30_masked_fit.rds (for the pseudobulk +
##   metadata) -> results/stephenson_K30_masked_fit_excl.rds

suppressPackageStartupMessages(library(cellprograms))

a <- commandArgs(TRUE)
fit_rds <- if (length(a) >= 1) a[1] else file.path("results", "stephenson_K30_masked_fit.rds")
out_rds <- if (length(a) >= 2) a[2] else file.path("results", "stephenson_K30_masked_fit_excl.rds")

f <- readRDS(fit_rds)
Y <- f$pseudobulk
meta <- f$sample_metadata
all_genes <- unique(unlist(lapply(Y, colnames)))
excl_old <- f$preprocessing$exclude_genes
excl_new <- artifact_genes(all_genes)                       # all 4 sets
cat("universe:", length(all_genes),
    "| old excl:", length(excl_old),
    "| new excl (all sets):", length(excl_new),
    "| delta:", length(setdiff(excl_new, excl_old)), "\n")

x <- as_cell_program_data(Y, sample_metadata = meta)
t0 <- Sys.time()
fit <- canonicalize_programs(fit_celltype_programs(
  x, loading_prior = "point_laplace", var_type = 1, features = "variable",
  n_variable_genes = 2000L, exclude_artifacts = TRUE, max_factors = 30L,
  covariates = c("institute", "sex"), covariate_mode = "residualize", seed = 1L))
cat("fit time:", format(Sys.time() - t0), "\n")
stopifnot(isTRUE(fit$preprocessing$exclude_artifacts),
          setequal(fit$preprocessing$exclude_genes, excl_new))
fit$pseudobulk <- Y
saveRDS(fit, out_rds)
cat("programs:", sum(vapply(fit$loadings, ncol, 1L)), "\n")
