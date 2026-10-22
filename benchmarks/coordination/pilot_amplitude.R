## pilot_amplitude.R — oracle detectability floor vs amplitude, per breadth.
## Oracle = projection of adjusted matrices onto TRUE loadings, then the frozen
## node_cor machinery (stratified permutation, max-T FWER, alpha=0.05).
## Decision rule: per breadth, smallest amp with oracle pairwise power >= ~0.5.
## Sharded: SHARD_I / SHARD_N env vars; writes out/amp_shard_<i>.csv.
## run: SHARD_I=0 SHARD_N=8 R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/pilot_amplitude.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R")) source(file.path(here, f))

G <- as.integer(Sys.getenv("PILOT_G", "2000"))
options(ms.G = G)
shard_i <- as.integer(Sys.getenv("SHARD_I", "0"))
shard_n <- as.integer(Sys.getenv("SHARD_N", "1"))
dir.create(file.path(here, "out"), showWarnings = FALSE)

amps <- as.numeric(strsplit(Sys.getenv("AMPS", "0.2,0.3,0.5,0.75,1.0,1.5"), ",")[[1]])
tag <- Sys.getenv("OUT_TAG", "amp")
breadths <- c("R100", "R50", "R20", "R20small")
nrep <- c(R100 = 8L, R50 = 8L, R20 = 24L, R20small = 24L)

grid <- do.call(rbind, lapply(breadths, function(b)
  data.frame(breadth = b, amp = amps, rep = rep(seq_len(nrep[[b]]), each = length(amps)),
             stringsAsFactors = FALSE)))
grid <- grid[order(grid$breadth, grid$amp, grid$rep), ]
grid$seed <- 91000L + match(grid$breadth, breadths) * 10000L +
             match(grid$amp, amps) * 1000L + grid$rep
sel <- which(seq_len(nrow(grid)) %% shard_n == shard_i)

rows <- list()
for (gi in sel) {
  g <- grid[gi, ]
  d <- gen_coord("act_only", amp = g$amp, seed = g$seed, breadth = g$breadth)
  mats <- adjusted_matrices(d$x, "institute")
  o <- method_oracle(mats, d$truth, d$strata, n_perm = 999L, seed = g$seed)
  ed <- o$edges
  for (p in d$truth$programs) {
    mem <- p$members
    prs <- if (length(mem) >= 2L) combn(mem, 2L, simplify = FALSE) else list()
    n_found <- 0L
    if (!is.null(ed) && nrow(ed)) {
      hit <- ed$prog_a == p$id & ed$prog_b == p$id &
             ed$ct_a %in% mem & ed$ct_b %in% mem
      n_found <- sum(hit)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      breadth = g$breadth, amp = g$amp, rep = g$rep, program = p$id,
      n_pairs = length(prs), n_found = n_found,
      n_spurious = if (is.null(ed)) 0L else sum(!(ed$prog_a == ed$prog_b &
                    ed$prog_a %in% names(d$truth$programs))),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("[shard %d] %s amp=%.2f rep=%d done\n", shard_i, g$breadth, g$amp, g$rep))
}
out <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(out, file.path(here, "out", sprintf("%s_shard_%02d.csv", tag, shard_i)),
          row.names = FALSE)
cat(sprintf("[shard %d] wrote %d rows\n", shard_i, nrow(out)))
