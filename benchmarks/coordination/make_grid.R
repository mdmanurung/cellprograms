## make_grid.R — build the bounded main-run job grid (out/main_grid.csv).
## Amplitudes come from out/amp_choice.csv (oracle pilot); replicate counts
## follow prereg-coord-v1 §8. Seeds are per-replicate and shared across arms
## (all arms run inside one replicate process).
## run: Rscript benchmarks/coordination/make_grid.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })

choice <- read.csv(file.path(here, "out", "amp_choice.csv"), stringsAsFactors = FALSE)
amp_of <- setNames(choice$amp, choice$breadth)

N_PRIMARY <- as.integer(Sys.getenv("N_PRIMARY", "32"))
N_ACTGENES <- as.integer(Sys.getenv("N_ACTGENES", N_PRIMARY))  # dev: act+genes tracks the primary count; confirmatory sets 0 (primary endpoint is act_only)
N_NULL <- as.integer(Sys.getenv("N_NULL", "16"))
N_PRIV <- as.integer(Sys.getenv("N_PRIV", "24"))
N_GENES <- as.integer(Sys.getenv("N_GENES", "16"))
N_R100 <- as.integer(Sys.getenv("N_R100", "8"))
N_SWEEP <- as.integer(Sys.getenv("N_SWEEP", "16"))

rows <- list(); k <- 0L
add <- function(mode, breadth, amp, act_cor, n) {
  for (r in seq_len(n)) {
    k <<- k + 1L
    rows[[k]] <<- data.frame(mode = mode, breadth = breadth, amp = amp,
                             act_cor = act_cor, rep = r,
                             seed = 1000L + k * 37L + r)
  }
}
## primary cells
for (b in c("R50", "R20")) {
  add("act_only",  b, amp_of[[b]], 1, N_PRIMARY)
  add("act+genes", b, amp_of[[b]], 1, N_ACTGENES)
}
## secondary
add("null",       "R50", amp_of[["R50"]], 1, N_NULL)   # amp unused by null
add("priv_corr",  "R20", amp_of[["R20"]], 1, N_PRIV)   # decoys ~ R20-like pairs
add("genes_only", "R50", amp_of[["R50"]], 1, N_GENES)
add("act_only",   "R100", amp_of[["R100"]], 1, N_R100) # sanity
add("act_only",   "R20small", amp_of[["R20small"]], 1, N_GENES) # secondary: small-ct pairs
## act_cor sweep on the primary family (act_only, each-alone)
for (ac in c(0.9, 0.7)) for (b in c("R50", "R20")) add("act_only", b, amp_of[[b]], ac, N_SWEEP)

grid <- do.call(rbind, rows)
write.csv(grid, file.path(here, "out", "main_grid.csv"), row.names = FALSE)
cat(sprintf("grid: %d replicates\n", nrow(grid)))
print(table(grid$mode, grid$breadth))
