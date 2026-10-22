## pilot_stage1.R — stage-1 recovery check at the oracle-chosen amplitudes.
## Reads out/amp_choice.csv (columns: breadth, amp). For each (breadth, amp,
## rep): stage-1 fit, then per (program, member ct) the max |cor| between a
## canonical program's scores and true activity. Recovery rule: |r| >= 0.5.
## Prereg rule: recovery fraction >= 0.5, else move one amp grid step up and
## record the deviation. Sharded like pilot_amplitude.R.
## run: SHARD_I=0 SHARD_N=4 R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/pilot_stage1.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R")) source(file.path(here, f))

G <- as.integer(Sys.getenv("PILOT_G", "2000"))
options(ms.G = G)
shard_i <- as.integer(Sys.getenv("SHARD_I", "0"))
shard_n <- as.integer(Sys.getenv("SHARD_N", "1"))
n_rep <- as.integer(Sys.getenv("STAGE1_REPS", "4"))
dir.create(file.path(here, "out"), showWarnings = FALSE)

choice <- read.csv(file.path(here, "out", "stage1_grid.csv"), stringsAsFactors = FALSE)
grid <- do.call(rbind, lapply(seq_len(nrow(choice)), function(i)
  data.frame(breadth = choice$breadth[i], amp = choice$amp[i], rep = seq_len(n_rep))))
grid$seed <- 95000L + match(grid$breadth, c("R100", "R50", "R20", "R20small")) * 1000L +
             round(grid$amp * 100) + grid$rep
sel <- which(seq_len(nrow(grid)) %% shard_n == shard_i)

rows <- list()
for (gi in sel) {
  g <- grid[gi, ]
  d <- gen_coord("act_only", amp = g$amp, seed = g$seed, breadth = g$breadth)
  fit <- fit_celltype_programs(d$x, max_factors = 10L, covariates = "institute",
                               covariate_mode = "residualize", seed = g$seed)
  fit <- canonicalize_programs(fit)
  for (p in d$truth$programs) {
    for (ct in p$members) {
      S <- fit$scores[[ct]]
      zt <- p$z[[ct]][rownames(S)]
      best <- if (ncol(S)) max(abs(cor(S, zt))) else 0
      rows[[length(rows) + 1L]] <- data.frame(
        breadth = g$breadth, amp = g$amp, rep = g$rep, program = p$id, cell_type = ct,
        n_programs_ct = ncol(S), max_abs_cor = best, recovered = best >= 0.5,
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("[shard %d] %s amp=%.2f rep=%d done\n", shard_i, g$breadth, g$amp, g$rep))
}
out <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(out, file.path(here, "out", sprintf("stage1_shard_%02d.csv", shard_i)),
          row.names = FALSE)
cat(sprintf("[shard %d] wrote %d rows\n", shard_i, nrow(out)))
