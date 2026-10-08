## pilot_amp.R -- apply the PREREG amplitude rule to results/pilot: base power closest to 0.5, inside [0.3, 0.7] if possible.
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
f <- list.files(file.path(here, "results", "pilot"), pattern = "_pair\\.csv$", recursive = TRUE, full.names = TRUE)
f <- f[!grepl("runtime_", f)]
d <- do.call(rbind, lapply(f, function(p) read.csv(p, stringsAsFactors = FALSE)))
d <- d[d$arm == "base" & d$is_true, ]
tab <- aggregate(shared ~ base + scenario + amp, d, function(v) c(P = mean(v), n = length(v)))
tab <- do.call(data.frame, tab); names(tab)[4:5] <- c("P", "n")
tab12 <- tab[tab$scenario %in% c("S1", "S2"), ]
pick <- do.call(rbind, lapply(split(tab12, list(tab12$base, tab12$scenario)), function(g) {
  ins <- g[g$P >= 0.3 & g$P <= 0.7, ]
  h <- if (nrow(ins)) ins else g
  h <- h[order(abs(h$P - 0.5), -h$amp), ]     # ties go to the larger amplitude
  h[1, c("base", "scenario", "amp", "P")]
}))
## DEVIATIONS D1: S3 base power is capped by callability, so S3 uses the S2 amplitude (neutral to arms); S7 uses S1.
pick <- pick[pick$scenario %in% c("S1", "S2"), ]
for (b in unique(pick$base)) pick <- rbind(pick, transform(pick[pick$base == b & pick$scenario == "S2", ], scenario = "S3", P = NA))
print(tab[order(tab$base, tab$scenario, tab$amp), ], row.names = FALSE); print(pick, row.names = FALSE)
write.csv(pick[, c("base", "scenario", "amp")], file.path(here, "amplitudes.csv"), row.names = FALSE)
