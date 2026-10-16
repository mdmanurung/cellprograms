## pilot.R -- pilot for prereg-coord-v1: (1) calibrate the amplitude per breadth class so the
## pairwise oracle (true-z correlation test on fitted scores) has power about 0.5, (2) pick the
## node density filter margin, (3) measure per-replicate fit cost for the main-run estimate.
##
## Sandbox smoke:  Rscript pilot.R --reps 2 --amps 0.4,0.8,1.6 --nperm 199
## Full pilot:     Rscript pilot.R --reps 10 --amps 0.2,0.4,0.8,1.6,3.2 --nperm 999
##
## Outputs outputs/pilot.csv (per amp x breadth: oracle power, per-replicate fit seconds) and
## outputs/pilot_density.csv (null-scenario cluster counts per min_density).

args <- commandArgs(TRUE)
get_arg <- function(flag, default) {
  i <- match(paste0("--", flag), args)
  if (is.na(i)) return(default) else args[i + 1L]
}
reps <- as.integer(get_arg("reps", "2"))
amps <- as.numeric(strsplit(get_arg("amps", "0.4,0.8,1.6"), ",")[[1]])
nperm <- as.integer(get_arg("nperm", "199"))
options(ms.G = as.integer(get_arg("G", "4000")))
dir.create("outputs", showWarnings = FALSE)

here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
source(file.path(here, "sim_coord.R"))
source(file.path(here, "methods.R"))
source(file.path(here, "scoring.R"))

FIT <- function(x, seed) {
  canonicalize_programs(fit_celltype_programs(
    x, loading_prior = "point_laplace", max_factors = 30, center = TRUE, scale = FALSE,
    features = "all", var_type = 1, backfit = TRUE, nullcheck = TRUE,
    covariates = ~institute, covariate_mode = "residualize", seed = seed))
}

## pairwise oracle on fitted scores: for a planted program, take the member cell-type pair with
## the highest |cor| between the nodes best-matched to the true z, and test it against the
## within-strata permutation null (the same statistic node_cor uses).
oracle_power <- function(fit, truth, strata, n_perm) {
  hits <- 0; tot <- 0
  for (s in names(truth$programs)) {
    p <- truth$programs[[s]]; m <- p$members
    if (length(m) < 2L) next
    best <- 0; best_va <- NULL; best_vb <- NULL
    for (i in seq_len(length(m) - 1L)) for (j in (i + 1L):length(m)) {
      a <- m[i]; b <- m[j]
      ids <- intersect(rownames(fit$scores[[a]]), rownames(fit$scores[[b]]))
      if (length(ids) < 20L) next
      za <- p$z[[a]][ids]; zb <- p$z[[b]][ids]
      va <- fit$scores[[a]][ids, which.max(abs(cor(fit$scores[[a]][ids, , drop = FALSE], za)))]
      vb <- fit$scores[[b]][ids, which.max(abs(cor(fit$scores[[b]][ids, , drop = FALSE], zb)))]
      r <- abs(suppressWarnings(cor(va, vb)))
      if (r > best) { best <- r; best_va <- va; best_vb <- vb }
    }
    if (best == 0) next
    set.seed(42)   # fixed permutation seed for the oracle across amps
    mx <- replicate(n_perm, abs(cor(best_va[sample.int(length(best_va))],
                                    best_vb[sample.int(length(best_vb))])))
    tot <- tot + 1L
    hits <- hits + as.integer((1 + sum(mx >= best)) / (n_perm + 1L) < 0.05)
  }
  c(hits = hits, tot = tot)
}

rows <- list()
for (amp in amps) {
  for (rep in seq_len(reps)) {
    seed <- 1000L * match(amp, amps) + rep
    d <- gen_coord("act+genes", amp = amp, seed = seed)
    t0 <- Sys.time()
    fit <- FIT(d$x, seed)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    op <- oracle_power(fit, d$truth, d$strata, nperm)
    rows[[length(rows) + 1L]] <- data.frame(amp = amp, rep = rep, fit_secs = secs,
                                            oracle_hits = op[["hits"]], oracle_tot = op[["tot"]])
    message(sprintf("amp %.1f rep %d: fit %.0f s, oracle %d/%d", amp, rep, secs, op[["hits"]], op[["tot"]]))
  }
}
pilot <- do.call(rbind, rows)
write.csv(pilot, "outputs/pilot.csv", row.names = FALSE)
pw <- aggregate(cbind(oracle_hits, fit_secs) ~ amp, pilot, sum)
pw$power <- pw$oracle_hits / reps
message("oracle power by amplitude (target ~0.5):")
message(paste(capture.output(print(pw[, c("amp", "power", "fit_secs")])), collapse = "\n"))

## density margin: null-scenario cluster counts under candidate min_density values
## (same seed per method x rep, so the edges are identical across densities and only the
## filter differs; three densities keep the permutation cost at 3x).
dens_rows <- list()
for (rep in seq_len(min(5L, reps))) {
  d <- gen_coord("null", amp = 1, seed = 500L + rep)
  fit <- FIT(d$x, 500L + rep)
  for (m in c("node_cor", "node_agg")) {
    for (md in c(0, 0.5, 0.9)) {
      r <- run_method(m, fit, d$strata, n_perm = nperm, seed = 1, min_density = md)
      dens_rows[[length(dens_rows) + 1L]] <- data.frame(method = m, min_density = md, rep = rep,
                                                        n_clusters = length(unique(r$clusters$cluster)))
    }
  }
  message(sprintf("null rep %d done", rep))
}
dens <- do.call(rbind, dens_rows)
write.csv(dens, "outputs/pilot_density.csv", row.names = FALSE)
message("null cluster counts by min_density (pick the smallest removing >= 90% of null clusters):")
message(paste(capture.output(aggregate(n_clusters ~ method + min_density, dens, mean)), collapse = "\n"))
