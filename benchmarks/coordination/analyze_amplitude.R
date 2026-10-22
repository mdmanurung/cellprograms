## analyze_amplitude.R — oracle power curves -> per-breadth amplitude choice.
## Rule (prereg §4): smallest amp with oracle pairwise power >= 0.5.
## Flags: "extend_down" if power at grid floor > 0.6 (curve not bracketed);
## "ceiling" if no amp reaches 0.5 (take max amp, record deviation).
## Writes out/amp_choice.csv and prints the full power table.
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/analyze_amplitude.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })

fs <- list.files(file.path(here, "out"),
                 pattern = sprintf("^%s_shard_.*\\.csv$", Sys.getenv("OUT_TAG", "amp")),
                 full.names = TRUE)
stopifnot("no amp shards found" = length(fs) > 0L)
dat <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
cat(sprintf("loaded %d rows from %d shards\n", nrow(dat), length(fs)))

pw <- aggregate(cbind(n_found, n_pairs) ~ breadth + amp, dat, sum)
pw$power <- pw$n_found / pw$n_pairs
pw$se <- sqrt(pw$power * (1 - pw$power) / pw$n_pairs)
spur <- aggregate(n_spurious ~ breadth + amp, dat, mean)
pw <- merge(pw, spur, by = c("breadth", "amp"))
pw <- pw[order(pw$breadth, pw$amp), ]
cat("\n== oracle pairwise power (act_only) ==\n")
print(pw, row.names = FALSE, digits = 2)

choice <- do.call(rbind, lapply(split(pw, pw$breadth), function(g) {
  g <- g[order(g$amp), ]
  hit <- which(g$power >= 0.5)
  flag <- ""
  if (!length(hit)) {
    pick <- nrow(g); flag <- "ceiling: max amp taken, power < 0.5 (deviation)"
  } else {
    pick <- hit[1]
    if (pick == 1L && g$power[1] > 0.6) flag <- "extend_down: grid floor already > 0.6"
  }
  data.frame(breadth = g$breadth[1], amp = g$amp[pick],
             oracle_power = g$power[pick], flag = flag, stringsAsFactors = FALSE)
}))
cat("\n== amplitude choice ==\n")
print(choice, row.names = FALSE)
write.csv(choice, file.path(here, "out", "amp_choice.csv"), row.names = FALSE)
