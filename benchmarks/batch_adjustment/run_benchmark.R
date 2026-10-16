#!/usr/bin/env Rscript
# batch-adj-v1 harness (PREREG.md in this directory). Runs the four adjustment
# arms on the batch-planted simulator and writes per-seed results RDS.
# Usage:
#   Rscript run_benchmark.R --sigma 0.375 --nseed 10 --out results_s0375 \
#           [--confound 0] [--stability 1]
library(cellprograms)

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1L]
}
sigma_b <- as.numeric(getarg("--sigma", 0.375))
nseed <- as.integer(getarg("--nseed", 10))
outdir <- getarg("--out", paste0("results_s", gsub("\\.", "", sigma_b)))
confound_a <- as.numeric(getarg("--confound", 0))
n_samples <- as.integer(getarg("--samples", 150))
do_stab <- as.integer(getarg("--stability", 1)) == 1L
cts <- unlist(strsplit(getarg("--cts", "B,CD8"), ","))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "sim_batch.R"))

# --- metrics -----------------------------------------------------------------
recovery <- function(fit, ct, truth) {
  S <- fit$scores[[ct]]
  if (is.null(S) || ncol(S) == 0L) return(NA_real_)
  zs <- truth[[ct]]
  mean(vapply(zs, function(p) {
    if (ncol(S) < 1L) return(NA_real_)
    max(abs(cor(S[, seq_len(ncol(S)), drop = FALSE], p$z)))
  }, numeric(1)))
}

batch_r2 <- function(fit, ct, batch) {
  S <- fit$scores[[ct]]
  if (is.null(S) || ncol(S) == 0L) return(NA_real_)
  mean(vapply(seq_len(ncol(S)), function(k) {
    summary(lm(S[, k] ~ batch))$r.sq
  }, numeric(1)))
}

# --- one seed, one cell type, four arms --------------------------------------
run_seed_ct <- function(seed, ct) {
  sim <- simulate_pseudobulk(scenario = "mixed", n_samples = n_samples,
                             n_genes = 1500, program_size = 50,
                             signal_sd = 1.5, noise_sd = 1.0, seed = seed)
  bt <- plant_batch(sim, n_batches = 4, sigma_b = sigma_b,
                    confound_a = confound_a, seed = seed + 100L)
  meta <- data.frame(observation_id = names(bt$batch), batch = bt$batch)
  rownames(meta) <- meta$observation_id
  out <- list(seed = seed, ct = ct, sigma_b = sigma_b, confound_a = confound_a,
              n_samples = n_samples)

  # A0: OLS residualization (default)
  x0 <- as_cell_program_data(bt$pseudobulk[ct], sample_metadata = meta)
  f0 <- tryCatch(canonicalize_programs(fit_celltype_programs(
    x0, max_factors = 10, covariates = ~batch, seed = 1)),
    error = function(e) structure(list(err = conditionMessage(e)), class = "fit_error"))
  # A3: within-EBMF fixed factors
  x3 <- as_cell_program_data(bt$pseudobulk[ct], sample_metadata = meta)
  f3 <- tryCatch(canonicalize_programs(fit_celltype_programs(
    x3, max_factors = 10, covariates = ~batch, covariate_mode = "fixed", seed = 1)),
    error = function(e) structure(list(err = conditionMessage(e)), class = "fit_error"))
  # A1: ComBat on log-pseudobulk, then no covariates
  Y <- bt$pseudobulk[[ct]]
  Yc <- t(sva::ComBat(t(Y), batch = bt$batch))
  dimnames(Yc) <- dimnames(Y)
  x1 <- as_cell_program_data(stats::setNames(list(Yc), ct))
  f1 <- tryCatch(canonicalize_programs(fit_celltype_programs(
    x1, max_factors = 10, seed = 1)),
    error = function(e) structure(list(err = conditionMessage(e)), class = "fit_error"))
  # A2: ComBat-seq on count-pseudobulk, re-normalized, then no covariates
  counts <- make_counts(bt$pseudobulk[ct], seed = seed + 200L)[[ct]]
  adj <- t(sva::ComBat_seq(counts = t(counts), batch = bt$batch))
  dimnames(adj) <- dimnames(counts)
  Y2 <- lognorm_counts(adj)
  x2 <- as_cell_program_data(stats::setNames(list(Y2), ct))
  f2 <- tryCatch(canonicalize_programs(fit_celltype_programs(
    x2, max_factors = 10, seed = 1)),
    error = function(e) structure(list(err = conditionMessage(e)), class = "fit_error"))

  arms <- list(A0 = f0, A1 = f1, A2 = f2, A3 = f3)
  for (nm in names(arms)) {
    f <- arms[[nm]]
    if (inherits(f, "fit_error")) {
      out[[nm]] <- list(error = f$err)
      next
    }
    rec <- recovery(f, ct, sim$truth)
    br2 <- batch_r2(f, ct, bt$batch)
    K <- ncol(f$scores[[ct]])
    stab <- NA_real_
    if (do_stab && identical(sigma_b, 0.375) && identical(ct, "B")) {
      s <- tryCatch(assess_program_stability(f, n_boot = 5L, seed = 1,
                                             progress = FALSE), error = function(e) NULL)
      if (!is.null(s) && !is.null(s$stability)) {
        stab <- mean(s$stability$per_program$median_score_cor, na.rm = TRUE)
      }
    }
    out[[nm]] <- list(recovery = rec, batch_r2 = br2, K = K, stability = stab)
  }
  out
}

results <- list()
for (seed in seq_len(nseed)) {
  for (ct in cts) {
    r <- tryCatch(run_seed_ct(seed, ct), error = function(e) list(
      seed = seed, ct = ct, fatal = conditionMessage(e)))
    results[[length(results) + 1L]] <- r
    saveRDS(results, file.path(outdir, "results_partial.rds"))
    cat(sprintf("seed %d ct %s done\n", seed, ct))
  }
}
saveRDS(results, file.path(outdir, "results.rds"))
cat("DONE", length(results), "runs ->", file.path(outdir, "results.rds"), "\n")
