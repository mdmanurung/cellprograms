## pilot_runtime.R — time every race arm once at candidate full scale.
## Decides G and replicate counts for the bounded main run (block_E is the
## expected bottleneck). Run EXCLUSIVELY (no other R processes) or timings lie.
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/pilot_runtime.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R")) source(file.path(here, f))

G <- as.integer(Sys.getenv("PILOT_G", "2000"))
options(ms.G = G)
cat(sprintf("== runtime pilot at G=%d ==\n", G))

d <- gen_coord("act_only", amp = 0.5, seed = 777, breadth = "R50")
tim <- list()

t0 <- Sys.time()
fit <- fit_celltype_programs(d$x, max_factors = 10L, covariates = "institute",
                             covariate_mode = "residualize", seed = 777)
fit <- canonicalize_programs(fit)
tim$stage1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("stage-1 (flashier, K=10): %.0fs | programs/ct: %s\n", tim$stage1,
            paste(vapply(fit$scores, ncol, 1L), collapse = ",")))

for (m in c("pa_cc", "node_cor", "subspace_sum", "gene_match", "block_E", "eb_cp")) {
  gc()
  t0 <- Sys.time()
  res <- run_method(m, fit, d$strata, seed = 1L)
  tim[[m]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("%-12s %7.0fs | %d calls\n", m, tim[[m]], length(unique(res$clusters$cluster))))
}

tot <- sum(unlist(tim))
cat(sprintf("\nTOTAL per dataset: %.0fs (%.1f min)\n", tot, tot / 60))
cat(sprintf("Extrapolation, 192 datasets, 12 parallel workers: %.1f h\n",
            tot * 192 / 12 / 3600))
saveRDS(list(G = G, timings = tim, total_secs = tot),
        file.path(here, "out", sprintf("pilot_runtime_G%d.rds", G)))
