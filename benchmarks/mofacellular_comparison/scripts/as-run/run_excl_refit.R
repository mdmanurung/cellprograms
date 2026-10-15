## run_excl_refit.R -- sensitivity refit: exclude_artifacts = TRUE (all 4 sets).
## Same recipe as the masked K30 refit (point_laplace, K=30, var_type=1, backfit,
## nullcheck, covariates institute+sex, residualize, 2000 HVGs, seed 1) but with
## the full artifact preset instead of the vdj+genotype_marker 274-gene list.
## Exclusion still happens before HVG selection, so fits select 2000 HVGs.
.libPaths("/workspace/.Rlib")
pkgload::load_all("/workspace/repo2", quiet = TRUE)

f <- readRDS("/workspace/stephenson_K30_masked_fit.rds")
Y <- f$pseudobulk
meta <- f$sample_metadata
all_genes <- unique(unlist(lapply(Y, colnames)))
excl_old <- f$preprocessing$exclude_genes
excl_new <- artifact_genes(all_genes)                       # all 4 sets
cat("universe:", length(all_genes),
    "| old excl (vdj+genotype_marker):", length(excl_old),
    "| new excl (all sets):", length(excl_new),
    "| delta:", length(setdiff(excl_new, excl_old)), "\n")
cat("delta genes:", paste(sort(setdiff(excl_new, excl_old)), collapse = ", "), "\n")
flush.console()

x <- as_cell_program_data(Y, sample_metadata = meta)
t0 <- Sys.time()
fit <- canonicalize_programs(fit_celltype_programs(
  x, loading_prior = "point_laplace", var_type = 1, features = "variable",
  n_variable_genes = 2000L, exclude_artifacts = TRUE, max_factors = 30L,
  covariates = c("institute", "sex"), covariate_mode = "residualize", seed = 1L))
cat("fit time:", format(Sys.time() - t0), "\n"); flush.console()
stopifnot(isTRUE(fit$preprocessing$exclude_artifacts),
          setequal(fit$preprocessing$exclude_genes, excl_new))
fit$pseudobulk <- Y
saveRDS(fit, "/workspace/stephenson_K30_masked_fit_excl.rds")
cat("programs:", sum(vapply(fit$loadings, ncol, 1L)), "\n"); flush.console()

t0 <- Sys.time()
cmp <- find_coordinated_programs(fit, n_perm = 999, seed = 1, min_donors = 20)
cat("CMP time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(cmp, "/workspace/stephenson_K30_cmp_excl.rds")
cat("n CMPs:", length(unique(cmp$cmps$cmp)), "| edges:", nrow(cmp$edges), "\n")
cat("ALL DONE\n")
