## paired comparison of two candidates on per-dataset mean target F1 (marginal mode, positive scenarios)
args <- commandArgs(TRUE); dir <- args[1]; a <- args[2]; b <- args[3]
f <- function(cn) { g <- as.data.frame(readRDS(file.path(dir, paste0(cn, ".rds")))$signals); g <- g[g$mode == "marginal", ]
  stats::aggregate(f1 ~ base + scenario + seed, g, mean) }
A <- f(a); B <- f(b); m <- merge(A, B, by = c("base", "scenario", "seed")); d <- m$f1.y - m$f1.x
set.seed(1); bs <- replicate(5000, mean(sample(d, replace = TRUE)))
cat(sprintf("%s - %s: mean diff %.4f  95%% CI [%.4f, %.4f]  n=%d  better/worse/tie %d/%d/%d  wilcox p=%.3g\n", b, a, mean(d),
  quantile(bs, .025), quantile(bs, .975), length(d), sum(d > 0), sum(d < 0), sum(d == 0), suppressWarnings(stats::wilcox.test(d)$p.value)))
for (bb in unique(m$base)) { i <- m$base == bb; cat("  ", bb, sprintf("%.4f", mean(d[i])), "\n") }
