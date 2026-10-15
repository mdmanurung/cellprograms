## run_sharing.R -- program_sharing() on the real Stephenson/COMBAT K30 fits.
## The real fits already store F_psd/L_psd (full flashier object kept), column-aligned with the
## canonicalized loadings (verified: F_pm * sign * scale == loadings exactly, rownames identical),
## so no refit is needed. Runs the shipped program_sharing() (permutation-calibrated matching,
## CorShrink-style V, mashr) per cohort, plus flag_artifact_genes() on every program for the
## artifact overlay. Saves per-cohort .rds + program-level flag CSV to /workspace.
.libPaths("/workspace/.Rlib")
pkgload::load_all("/workspace/repo2", quiet = TRUE)
stopifnot(requireNamespace("mashr", quietly = TRUE))
R <- "/mnt/shared-workspace/cmp/real"
for (co in c("stephenson", "combat")) {
  fit <- readRDS(file.path(R, sprintf("%s_pb_sex_K30_fit.rds", co)))
  np <- sum(vapply(fit$loadings, ncol, 1L))
  cat("==", co, "| cell types:", length(fit$cell_types), "| programs:", np, "\n"); flush.console()
  t0 <- Sys.time()
  res <- program_sharing(fit, n_perm = 199, seed = 1, z_thresh = 1.5)
  cat("sharing runtime:", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n"); flush.console()
  saveRDS(res, sprintf("/workspace/%s_K30_sharing.rds", co))
  fl <- flag_artifact_genes(fit)
  write.csv(fl, sprintf("/workspace/%s_K30_prog_artifact_flags.csv", co), row.names = FALSE)
  cat("artifact-flagged programs:", sum(apply(fl[, grepl("_flag$", names(fl)), drop = FALSE], 1, any)), "/", nrow(fl), "\n")
  print(res); flush.console()
}
cat("ALL DONE\n")
