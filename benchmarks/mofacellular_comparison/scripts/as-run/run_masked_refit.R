## run_masked_refit.R -- B1: artifact-masked K30 refits + sharing rerun.
## Same recipe as the unmasked real fits (point_laplace, K=30, var_type=1, backfit,
## nullcheck, covariates institute+sex, residualize, 2000 HVGs, seed 1) plus
## exclude_genes = artifact_genes(all genes, sets = c("vdj", "genotype_marker")).
## Exclusion happens before HVG selection, so fits still select 2000 HVGs.
## Then program_sharing(n_perm=199, seed=1, z_thresh=1.5) -- identical settings to
## the unmasked run for comparability. Checkpoints each cohort to
## /mnt/shared-workspace/cmp/real/ immediately after completion.
.libPaths("/workspace/.Rlib")
pkgload::load_all("/workspace/repo2", quiet = TRUE)
stopifnot(requireNamespace("mashr", quietly = TRUE))
R <- "/mnt/shared-workspace/cmp/real"

for (co in c("stephenson", "combat")) {
  cat("==", co, "==\n"); flush.console()
  fit0 <- readRDS(file.path(R, sprintf("%s_pb_sex_K30_fit.rds", co)))
  Y <- fit0$pseudobulk
  meta <- fit0$sample_metadata
  all_genes <- unique(unlist(lapply(Y, colnames)))
  excl <- artifact_genes(all_genes, sets = c("vdj", "genotype_marker"))
  cat("gene universe:", length(all_genes), "| excluded:", length(excl), "\n"); flush.console()

  x <- as_cell_program_data(Y, sample_metadata = meta)
  t0 <- Sys.time()
  fit <- canonicalize_programs(fit_celltype_programs(
    x, loading_prior = "point_laplace", var_type = 1, features = "variable",
    n_variable_genes = 2000L, exclude_genes = excl, max_factors = 30L,
    covariates = c("institute", "sex"), covariate_mode = "residualize", seed = 1L))
  cat("fit time:", format(Sys.time() - t0), "\n"); flush.console()
  fit$pseudobulk <- Y
  saveRDS(fit, sprintf("/workspace/%s_K30_masked_fit.rds", co))
  system2("cp", c(sprintf("/workspace/%s_K30_masked_fit.rds", co),
                  file.path(R, sprintf("%s_K30_masked_fit.rds", co))))
  cat("programs:", sum(vapply(fit$loadings, ncol, 1L)), "\n"); flush.console()

  t0 <- Sys.time()
  res <- program_sharing(fit, n_perm = 199, seed = 1, z_thresh = 1.5)
  cat("sharing runtime:", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n"); flush.console()
  saveRDS(res, sprintf("/workspace/%s_K30_masked_sharing.rds", co))
  system2("cp", c(sprintf("/workspace/%s_K30_masked_sharing.rds", co),
                  file.path(R, sprintf("%s_K30_masked_sharing.rds", co))))

  fl <- flag_artifact_genes(fit)
  write.csv(fl, sprintf("/workspace/%s_K30_masked_prog_artifact_flags.csv", co), row.names = FALSE)
  flagcols <- grep("_flag$", names(fl), value = TRUE)
  cat("artifact-flagged programs (masked fit):", sum(apply(fl[, flagcols, drop = FALSE], 1, any)),
      "/", nrow(fl), "\n")
  print(res); flush.console()
}
cat("ALL DONE\n")
