## choose_amplitude.R — final per-breadth amplitude for prereg-coord-v1.
## Rule (prereg §4): oracle pick = smallest amp with oracle pairwise power >= 0.5
## (both oracle grids merged). Stage-1 veto: recovery = fraction of (program,
## member-ct) cells with max |cor(scores, true activity)| >= 0.5 must be >= 0.5;
## final amp = max(oracle pick, stage-1 minimum). Deviations recorded in the
## output table. Writes out/amp_choice.csv.
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/choose_amplitude.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })

read_pat <- function(pat) {
  fs <- list.files(file.path(here, "out"), pattern = pat, full.names = TRUE)
  if (!length(fs)) return(data.frame())
  do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
}
orc <- rbind(read_pat("^amp_shard_.*\\.csv$"), read_pat("^amplow_shard_.*\\.csv$"))
s1  <- read_pat("^stage1_shard_.*\\.csv$")
stopifnot(nrow(orc) > 0, nrow(s1) > 0)

pw <- aggregate(cbind(n_found, n_pairs) ~ breadth + amp, orc, sum)
pw$oracle_power <- pw$n_found / pw$n_pairs
rec <- aggregate(recovered ~ breadth + amp, s1, mean)
m <- merge(pw, rec, by = c("breadth", "amp"), all = TRUE)
m <- m[order(m$breadth, m$amp), ]
cat("== combined pilot table ==\n")
print(m[, c("breadth", "amp", "oracle_power", "recovered")], row.names = FALSE, digits = 2)

choice <- do.call(rbind, lapply(split(m, m$breadth), function(g) {
  g <- g[order(g$amp), ]
  o <- g$amp[which(g$oracle_power >= 0.5)[1]]
  s <- g$amp[which(g$recovered >= 0.5)[1]]
  final <- max(o, s, na.rm = TRUE)
  data.frame(breadth = g$breadth[1], amp = final,
             oracle_pick = o, stage1_min = s,
             oracle_power_at_choice = g$oracle_power[g$amp == final],
             stage1_recovery_at_choice = g$recovered[g$amp == final],
             deviation = ifelse(is.na(o) || final > o,
                                "stage-1 minimum above oracle pick (estimation-bound regime)",
                                ""),
             stringsAsFactors = FALSE)
}))
cat("\n== final amplitude choice ==\n")
print(choice, row.names = FALSE)
write.csv(choice, file.path(here, "out", "amp_choice.csv"), row.names = FALSE)
write.csv(m, file.path(here, "out", "pilot_combined.csv"), row.names = FALSE)
