## main_run.R — bounded main run for prereg-coord-v1: shard driver.
## Loops over a job grid CSV (make_grid.R); replicate logic lives in
## replicate.R. Per-replicate CSVs make the run resumable and isolate failures.
## run: GRID=... SHARD_I=0 SHARD_N=12 R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/main_run.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R", "replicate.R")) source(file.path(here, f))

grid_path <- Sys.getenv("GRID")
stopifnot("GRID env var required" = nzchar(grid_path))
shard_i <- as.integer(Sys.getenv("SHARD_I", "0"))
shard_n <- as.integer(Sys.getenv("SHARD_N", "1"))
G <- as.integer(Sys.getenv("MAIN_G", "2000"))
outdir <- Sys.getenv("OUTDIR", file.path(here, "out", "main"))
dir.create(file.path(outdir, "prog"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "ds"), recursive = TRUE, showWarnings = FALSE)
grid <- read.csv(grid_path, stringsAsFactors = FALSE)
mv <- Sys.getenv("METHODS", "ALL")
methods <- if (identical(mv, "ALL")) METHODS_ALL else strsplit(mv, "[, ]+")[[1]]
sel <- which(seq_len(nrow(grid)) %% shard_n == shard_i)
set.seed(shard_i * 7919L + 13L)          # desynchronize block_E memory peaks across shards
sel <- sample(sel)
cat(sprintf("[shard %d/%d] %d jobs\n", shard_i, shard_n, length(sel)))
for (gi in sel) {
  g <- grid[gi, ]
  st <- tryCatch(
    run_one_replicate(g$mode, g$breadth, g$amp, g$rep, g$seed, g$act_cor, G, outdir,
                      methods = methods),
    error = function(e) sprintf("FAILED: %s", conditionMessage(e)))
  cat(sprintf("[shard %d] %s r%03d: %s\n", shard_i,
              sprintf("%s_%s_a%g_c%g", g$mode, g$breadth, g$amp, g$act_cor), g$rep, st))
}
