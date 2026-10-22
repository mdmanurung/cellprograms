## highk_run.R — P3c secondary scenario: does the pa_cc default path degrade
## when the private-program load per cell type is high (~25 privates x 10 genes
## per cell type, on top of the R100+R50+R20 recurring programs)?  This
## scenario was NOT part of the 272-dataset preregistered grid; results are
## exploratory scope notes, not confirmatory.  Deviations from the frozen main
## run (documented): G = 4000, stage-1 max_factors = 30, private_per_ct = 25,
## n_priv_genes = 10.  Arms: pa_cc, node_cor, subspace_sum, gene_match, oracle
## (block_E / eb_cp skipped: their frozen caps, Kmax 20 / 15, are structurally
## misspecified for ~28 true factors per cell type).  Realized stage-1 K per
## cell type is recorded per replicate.
## run: SHARD_I=0 SHARD_N=12 R_LIBS_USER=/workspace/.Rlib \
##      OMP_NUM_THREADS=1 OPENBLAST_NUM_THREADS=1 Rscript benchmarks/coordination/highk_run.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "score.R")) source(file.path(here, f))

METHODS_HIGHK <- c("pa_cc", "node_cor", "subspace_sum", "gene_match", "oracle")
G_HIGHK <- 4000L
MAX_F_HIGHK <- 30L

make_highk_grid <- function() {
  ## 16 signal (act_only, all breadths planted, amp 0.5) + 8 null; seeds in a
  ## range disjoint from the frozen main grid (1000-3000).
  g <- rbind(
    data.frame(mode = "act_only", breadth = "all", amp = 0.5, act_cor = 1,
               rep = seq_len(16L), seed = 9000L + (seq_len(16L) - 1L) * 38L),
    data.frame(mode = "null", breadth = "all", amp = 0.5, act_cor = 1,
               rep = seq_len(8L), seed = 9500L + (seq_len(8L) - 1L) * 38L))
  g$cell_id <- sprintf("%s_%s_a%g_c%g", g$mode, g$breadth, g$amp, g$act_cor)
  g
}

run_one_highk <- function(mode, breadth, amp, rep, seed, act_cor = 1, outdir) {
  options(ms.G = G_HIGHK)
  cell_id <- sprintf("%s_%s_a%g_c%g", mode, breadth, amp, act_cor)
  fprog <- file.path(outdir, "prog", sprintf("%s_r%03d.csv", cell_id, rep))
  fds   <- file.path(outdir, "ds",   sprintf("%s_r%03d.csv", cell_id, rep))
  fk    <- file.path(outdir, "k",    sprintf("%s_r%03d.csv", cell_id, rep))
  if (file.exists(fds)) return(invisible("skip"))
  d <- gen_coord(mode, amp = amp, seed = seed, act_cor = act_cor, breadth = breadth,
                 G = G_HIGHK, private_per_ct = 25L, n_priv_genes = 10L)

  t0 <- Sys.time()
  fit <- fit_celltype_programs(d$x, max_factors = MAX_F_HIGHK, covariates = "institute",
                               covariate_mode = "residualize", seed = seed)
  fit <- canonicalize_programs(fit)
  t_stage1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ## realized stage-1 K per cell type (diagnostic: does max_factors bind?)
  kdf <- data.frame(cell_id = cell_id, mode = mode, rep = rep, seed = seed,
                    cell_type = fit$cell_types,
                    n_factors = vapply(fit$fits[fit$cell_types],
                                       function(f) as.integer(f$flash$n_factors), integer(1)),
                    max_factors = MAX_F_HIGHK, stringsAsFactors = FALSE)
  write.csv(kdf, fk, row.names = FALSE)

  prog_rows <- list(); ds_rows <- list()
  for (m in METHODS_HIGHK) {
    t0 <- Sys.time()
    res <- tryCatch(
      if (m == "oracle") {
        method_oracle(lapply(fit$fits[fit$cell_types], function(f) f$Y_used),
                      d$truth, d$strata, n_perm = 999L, seed = seed)
      } else {
        run_method(m, fit, d$strata, seed = seed)
      }, error = function(e) e)
    t_m <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
      ds_rows[[m]] <- data.frame(cell_id = cell_id, mode = mode, breadth = breadth,
                                 amp = amp, act_cor = act_cor, rep = rep, seed = seed,
                                 method = m, error = conditionMessage(res),
                                 runtime_secs = t_m, stage1_secs = t_stage1,
                                 n_calls = NA_integer_, n_spurious = NA_integer_,
                                 n_decoy = NA_integer_, n_trap = NA_integer_,
                                 spurious_per_call = NA_real_, stringsAsFactors = FALSE)
      next
    }
    s <- score_method(res, d$truth)
    dl <- cbind(data.frame(cell_id = cell_id, mode = mode, breadth = breadth,
                           amp = amp, act_cor = act_cor, rep = rep, seed = seed,
                           error = NA_character_, runtime_secs = t_m,
                           stage1_secs = t_stage1, stringsAsFactors = FALSE),
                s$dataset_level)
    ds_rows[[m]] <- dl
    if (nrow(s$program_level)) {
      prog_rows[[m]] <- cbind(dl[, c("cell_id", "mode", "breadth", "amp",
                                     "act_cor", "rep", "seed", "method")],
                              s$program_level)
    }
  }
  prog_df <- if (length(prog_rows)) do.call(rbind, prog_rows) else data.frame(
    cell_id = character(), mode = character(), breadth = character(),
    amp = numeric(), act_cor = numeric(), rep = integer(), seed = integer(),
    method = character(), program_id = character(), n_true = integer(),
    matched_cluster = integer(), recall = numeric(), precision = numeric(),
    breadth_error = integer(), found = logical(), stringsAsFactors = FALSE)
  write.csv(prog_df, fprog, row.names = FALSE)
  write.csv(do.call(rbind, ds_rows), fds, row.names = FALSE)
  invisible("done")
}

shard_i <- as.integer(Sys.getenv("SHARD_I", "0"))
shard_n <- as.integer(Sys.getenv("SHARD_N", "1"))
outdir <- file.path(here, "out", "highk")
for (sub in c("prog", "ds", "k")) dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
grid <- make_highk_grid()
write.csv(grid, file.path(here, "out", "highk_grid.csv"), row.names = FALSE)
sel <- which(seq_len(nrow(grid)) %% shard_n == shard_i)
set.seed(shard_i * 7919L + 13L)
sel <- sample(sel)
cat(sprintf("[shard %d/%d] %d jobs\n", shard_i, shard_n, length(sel)))
for (gi in sel) {
  g <- grid[gi, ]
  st <- tryCatch(
    run_one_highk(g$mode, g$breadth, g$amp, g$rep, g$seed, g$act_cor, outdir),
    error = function(e) sprintf("FAILED: %s", conditionMessage(e)))
  cat(sprintf("[shard %d] %s r%03d: %s\n", shard_i, g$cell_id, g$rep, st))
}
