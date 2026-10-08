## summarize_sim.R -- mean +/- sd over seeds by scenario x missing_frac x arm x metric
## Usage: Rscript summarize_sim.R [indir=results/sim_v2] [outfile=results/pa_vs_pmd_sim_summary_v2.csv]
## Writes the long summary CSV; prints the AUROC table first (scenarios with both shared
## and unshared pairs: S3, S9), then TPR/FPR for all scenarios, then failures.
options(width = 200)
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
a <- commandArgs(trailingOnly = TRUE)
indir <- file.path(here, if (length(a) >= 1) a[1] else "results/sim_v2")
outf  <- file.path(here, if (length(a) >= 2) a[2] else "results/pa_vs_pmd_sim_summary_v2.csv")
files <- list.files(indir, pattern = "^seed.*\\.csv$", full.names = TRUE)
d <- do.call(rbind, lapply(files, read.csv))
key <- d[c("scenario", "missing_frac", "arm", "metric")]
agg <- function(f) aggregate(d["value"], key, f)
s <- merge(merge(agg(function(v) mean(v, na.rm = TRUE)), agg(function(v) sd(v, na.rm = TRUE)),
                 by = names(key), suffixes = c("_mean", "_sd")),
           agg(function(v) sum(!is.na(v))), by = names(key))
names(s)[ncol(s)] <- "n_seeds"
s <- s[order(s$scenario, s$missing_frac, s$arm, s$metric), ]
s$summary <- sprintf("%.3f +/- %.3f", s$value_mean, s$value_sd)
write.csv(s, outf, row.names = FALSE)

wide <- function(metrics) {
  x <- s[s$metric %in% metrics, ]
  w <- reshape(x[c("scenario", "missing_frac", "arm", "metric", "value_mean")],
               idvar = c("scenario", "missing_frac", "arm"), timevar = "metric", direction = "wide")
  names(w) <- sub("value_mean.", "", names(w), fixed = TRUE)
  w[order(w$scenario, w$missing_frac, w$arm), c("scenario", "missing_frac", "arm", intersect(metrics, names(w)))]
}
pr <- function(w) print(format(w, digits = 3), row.names = FALSE)
cat(sprintf("seeds: %d files in %s\n\n== AUROC / TPR@FPR<=0.05 (S3, S9: mean over seeds) ==\n", length(files), indir))
pr(wide(c("auroc_cont", "auroc_sig", "tpr05_cont", "tpr05_sig", "tpr", "fpr")) |>
     (\(w) w[!is.na(w$auroc_cont) | !is.na(w$auroc_sig), ])())
cat("\n== TPR / FPR / n_false_shared / failed / nuis_est_cor (all scenarios) ==\n")
pr(wide(c("tpr", "fpr", "n_false_shared", "failed", "nuis_est_cor")))
