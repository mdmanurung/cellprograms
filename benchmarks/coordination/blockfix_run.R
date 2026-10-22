## blockfix_run.R — bounded rerun of the block_E arm with permutation-calibrated
## block-mass filtering (post-prereg development iteration, decision-table
## follow-up). 104 datasets from out/blockfix_grid.csv (identical seeds to the
## frozen main grid), block_E_cal arm only. Joint-fit seed matches the main
## run's block_E arm (seed = dataset seed), so block_E_cal differs from the
## main run's raw block_E ONLY by the max-T permutation layer (n_perm = 199,
## alpha = 0.05, within-strata donor permutation per cell-type block).
## Same per-replicate CSV contract as replicate.R -> resumable, isolates
## failures. Driver mirrors main_run.R sharding.
## run: GRID=out/blockfix_grid.csv SHARD_I=0 SHARD_N=12 R_LIBS_USER=/workspace/.Rlib \
##      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript benchmarks/coordination/blockfix_run.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
for (f in c("sim_coord.R", "methods.R", "score.R")) source(file.path(here, f))

METHOD <- "block_E_cal"
N_PERM <- 199L
ALPHA  <- 0.05

run_one_blockfix <- function(mode, breadth, amp, rep, seed, act_cor = 1,
                             G = 2000L, outdir) {
  options(ms.G = as.integer(G))
  cell_id <- sprintf("%s_%s_a%g_c%g", mode, breadth, amp, act_cor)
  fprog <- file.path(outdir, "prog", sprintf("%s_r%03d.csv", cell_id, rep))
  fds   <- file.path(outdir, "ds",   sprintf("%s_r%03d.csv", cell_id, rep))
  if (file.exists(fds)) return(invisible("skip"))
  d <- gen_coord(mode, amp = amp, seed = seed, act_cor = act_cor, breadth = breadth)

  t0 <- Sys.time()
  fit <- fit_celltype_programs(d$x, max_factors = 10L, covariates = "institute",
                               covariate_mode = "residualize", seed = seed)
  fit <- canonicalize_programs(fit)
  t_stage1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0 <- Sys.time()
  res <- tryCatch(
    method_block_e(fit, strata = d$strata, n_perm = N_PERM, seed = seed,
                   alpha = ALPHA),
    error = function(e) e)
  t_m <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "error")) {
    dl <- data.frame(cell_id = cell_id, mode = mode, breadth = breadth,
                     amp = amp, act_cor = act_cor, rep = rep, seed = seed,
                     method = METHOD, error = conditionMessage(res),
                     runtime_secs = t_m, stage1_secs = t_stage1,
                     n_calls = NA_integer_, n_spurious = NA_integer_,
                     n_decoy = NA_integer_, n_trap = NA_integer_,
                     spurious_per_call = NA_real_, stringsAsFactors = FALSE)
    write.csv(dl, fds, row.names = FALSE)
    write.csv(dl[0, ], fprog, row.names = FALSE)
    return(invisible("error"))
  }
  res$method <- METHOD                       # relabel for the results CSVs
  s <- score_method(res, d$truth)
  dl <- cbind(data.frame(cell_id = cell_id, mode = mode, breadth = breadth,
                         amp = amp, act_cor = act_cor, rep = rep, seed = seed,
                         error = NA_character_, runtime_secs = t_m,
                         stage1_secs = t_stage1, stringsAsFactors = FALSE),
              s$dataset_level)
  write.csv(dl, fds, row.names = FALSE)
  pl <- if (nrow(s$program_level)) {
    cbind(dl[, c("cell_id", "mode", "breadth", "amp", "act_cor", "rep",
                 "seed", "method")], s$program_level)
  } else {
    dl[0, c("cell_id", "mode", "breadth", "amp", "act_cor", "rep", "seed",
            "method")]
  }
  write.csv(pl, fprog, row.names = FALSE)
  invisible("done")
}

grid_path <- Sys.getenv("GRID", file.path(here, "out", "blockfix_grid.csv"))
stopifnot("GRID grid file not found" = file.exists(grid_path))
shard_i <- as.integer(Sys.getenv("SHARD_I", "0"))
shard_n <- as.integer(Sys.getenv("SHARD_N", "1"))
outdir <- file.path(here, "out", "blockfix")
dir.create(file.path(outdir, "prog"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "ds"), recursive = TRUE, showWarnings = FALSE)
grid <- read.csv(grid_path, stringsAsFactors = FALSE)
sel <- which(seq_len(nrow(grid)) %% shard_n == shard_i)
set.seed(shard_i * 7919L + 13L)
sel <- sample(sel)
cat(sprintf("[shard %d/%d] %d jobs\n", shard_i, shard_n, length(sel)))
for (gi in sel) {
  g <- grid[gi, ]
  st <- tryCatch(
    run_one_blockfix(g$mode, g$breadth, g$amp, g$rep, g$seed, g$act_cor, 2000L, outdir),
    error = function(e) sprintf("FAILED: %s", conditionMessage(e)))
  cat(sprintf("[shard %d] %s r%03d: %s\n", shard_i,
              sprintf("%s_%s_a%g_c%g", g$mode, g$breadth, g$amp, g$act_cor), g$rep, st))
}
