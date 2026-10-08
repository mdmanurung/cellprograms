## summarize_eval.R -- cp_L0 vs baselines verdict table. One command: Rscript summarize_eval.R
## Reads results/eval/<cohort>_<config>_corrected.csv, writes results/eval/verdicts.csv, prints a compact report.
## Method (per cohort x config x covariate; covariate roles in verdicts_cohort_spec.R):
##  * metric: categorical -> variant `availcase`, null = floor_masked (falls back to `floor` if floor_masked is
##    NA/absent; the `null_used` column says which); continuous -> `biomni_spearman`, null = `floor`
##    (dimension-matched). `cv_ridge` rows (continuous) are reported in cv_* columns only, never in the verdict.
##  * excess = value - null. Baselines = composition, clr_composition, global_pca, perct_pca, grouped_pca
##    (`random` is a control, not a baseline). Best baseline is chosen by EXCESS, not raw value.
##  * verdict_best (stricter; = final `label`) compares cp_L0 excess to the best baseline excess; verdict_median
##    to the median baseline excess. Tie band +-TIE (0.05 > ~0.033 reconstruction noise of this benchmark).
##  * technical covariates (Panel B, lower is better) get label cp_cleaner / tie / cp_worse; best = MOST confounded
##    baseline (max excess). Excluded (degenerate/duplicate) and unclassified covariates are never tallied.
.script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}
D <- .script_dir()
source(file.path(D, "verdicts_cohort_spec.R"))
TIE <- 0.05
BASE <- c("composition", "clr_composition", "global_pca", "perct_pca", "grouped_pca")
CP <- "cp_L0"

files <- sort(Sys.glob(file.path(D, "results/eval/*_corrected.csv")))
stopifnot(length(files) > 0)

role <- function(cohort, cv) {  # -> list(group, class, reason)
  s <- COHORT_SPEC[[cohort]]
  if (is.null(s)) return(list(group = "unclassified", class = NA, reason = "cohort not in spec"))
  if (cv %in% names(s$excl)) return(list(group = "excluded", class = NA, reason = unname(s$excl[cv])))
  if (cv %in% s$tech) return(list(group = "technical", class = NA, reason = ""))
  cls <- c(bio_expr = "expression", bio_comp = "compositional", bio_other = "other")
  for (nm in names(cls)) if (cv %in% s[[nm]]) return(list(group = "biology", class = cls[[nm]], reason = ""))
  warning("unclassified covariate: ", cohort, "/", cv); list(group = "unclassified", class = NA, reason = "not in spec")
}
verdict <- function(delta, tech = FALSE) {  # delta = cp excess - baseline excess
  if (is.na(delta)) return(NA_character_)
  if (tech) c("cp_cleaner", "tie", "cp_worse")[1 + (abs(delta) <= TIE) + 2 * (delta > TIE)]
  else c("baseline", "tie", "cp")[1 + (abs(delta) <= TIE) + 2 * (delta > TIE)]
}

rows <- list()
for (f in files) {
  m <- regmatches(basename(f), regexec("^([a-z0-9]+)_(.+)_corrected\\.csv$", basename(f)))[[1]]
  cohort <- m[2]; config <- m[3]
  x <- read.csv(f, stringsAsFactors = FALSE)
  if (!"floor_masked" %in% names(x)) x$floor_masked <- NA_real_
  for (cv in unique(x$covariate)) {
    xc <- x[x$covariate == cv, ]
    typ <- xc$type[1]
    mv <- if (typ == "categorical") "availcase" else "biomni_spearman"
    g <- xc[xc$metric_variant == mv, ]
    g <- g[!duplicated(g$representation), ]  # robust to re-run appends
    if (typ == "categorical") {
      use_masked <- !is.na(g$floor_masked)
      g$null <- ifelse(use_masked, g$floor_masked, g$floor)
      null_used <- if (all(use_masked)) "floor_masked" else if (!any(use_masked)) "floor(fallback)" else "mixed"
    } else { g$null <- g$floor; null_used <- "floor" }
    g$excess <- g$value - g$null
    cp <- g[g$representation == CP, ]; bl <- g[g$representation %in% BASE, ]
    if (nrow(cp) != 1 || !nrow(bl)) { warning("skipping ", cohort, "/", config, "/", cv, ": no cp or baselines"); next }
    r <- role(cohort, cv); tech <- r$group == "technical"
    bi <- which.max(bl$excess)  # best (biology) / most confounded (technical)
    med <- stats::median(bl$excess)
    d_best <- cp$excess - bl$excess[bi]; d_med <- cp$excess - med
    cvx <- xc[xc$metric_variant == "cv_ridge", ]; cvx <- cvx[!duplicated(cvx$representation), ]
    cvx$excess <- cvx$value - cvx$floor
    cvcp <- cvx[cvx$representation == CP, ]; cvbl <- cvx[cvx$representation %in% BASE, ]
    rows[[length(rows) + 1]] <- data.frame(
      cohort, config, covariate = cv, group = r$group, class = r$class, reason = r$reason,
      type = typ, metric = mv, null_used,
      cp_value = cp$value, cp_null = cp$null, cp_excess = cp$excess,
      best_baseline = bl$representation[bi], best_excess = bl$excess[bi], median_excess = med,
      d_best = d_best, d_median = d_med,
      verdict_best = verdict(d_best, tech), verdict_median = verdict(d_med, tech),
      label = if (r$group %in% c("biology", "technical")) verdict(d_best, tech) else r$group,
      cv_cp_excess = if (nrow(cvcp)) cvcp$excess else NA_real_,
      cv_best_excess = if (nrow(cvbl)) max(cvbl$excess) else NA_real_,
      cv_median_excess = if (nrow(cvbl)) stats::median(cvbl$excess) else NA_real_,
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
out <- out[order(out$cohort, out$config, factor(out$group, c("biology", "technical", "excluded", "unclassified")), out$covariate), ]
write.csv(out, file.path(D, "results/eval/verdicts.csv"), row.names = FALSE)

## ------------------------------------------------------------------ report
f3 <- function(v) formatC(v, format = "f", digits = 3)
show <- function(d, cols) { d[cols] <- lapply(d[cols], function(v) if (is.numeric(v)) f3(v) else v); print(d, row.names = FALSE, right = FALSE) }
cat(sprintf("\nverdicts.csv: %d rows | tie band +-%.2f | label = verdict vs BEST baseline (by excess)\n", nrow(out), TIE))

bio <- out[out$group == "biology", ]
cat("\n== TALLY: biology covariates, cohort x config x class (cp / tie / baseline; vs best | vs median) ==\n")
tal <- do.call(rbind, lapply(split(bio, list(bio$cohort, bio$config, bio$class), drop = TRUE), function(g) {
  n <- function(v, l) sum(v == l)
  data.frame(cohort = g$cohort[1], config = g$config[1], class = g$class[1], n = nrow(g),
             best = sprintf("%d/%d/%d", n(g$verdict_best, "cp"), n(g$verdict_best, "tie"), n(g$verdict_best, "baseline")),
             median = sprintf("%d/%d/%d", n(g$verdict_median, "cp"), n(g$verdict_median, "tie"), n(g$verdict_median, "baseline")))
}))
print(tal, row.names = FALSE, right = FALSE)
tot <- function(v) sprintf("%d/%d/%d", sum(v == "cp"), sum(v == "tie"), sum(v == "baseline"))
cat("\nOverall by config x class (best | median):\n")
print(do.call(rbind, lapply(split(bio, list(bio$config, bio$class), drop = TRUE), function(g)
  data.frame(config = g$config[1], class = g$class[1], n = nrow(g), best = tot(g$verdict_best), median = tot(g$verdict_median)))),
  row.names = FALSE, right = FALSE)

cat("\n== BIOLOGY covariates ==\n")
show(bio[c("cohort", "config", "covariate", "class", "cp_excess", "best_baseline", "best_excess",
           "median_excess", "d_best", "label", "verdict_median")],
     c("cp_excess", "best_excess", "median_excess", "d_best"))

cat("\n== TECHNICAL covariates, Panel B (cp excess vs MOST confounded baseline; lower is better) ==\n")
show(out[out$group == "technical", c("cohort", "config", "covariate", "cp_excess", "best_baseline",
                                      "best_excess", "median_excess", "d_best", "label", "verdict_median")],
     c("cp_excess", "best_excess", "median_excess", "d_best"))

cat("\n== EXCLUDED / UNCLASSIFIED (not tallied) ==\n")
ex <- unique(out[out$group %in% c("excluded", "unclassified"), c("cohort", "covariate", "reason")])
if (nrow(ex)) print(ex, row.names = FALSE, right = FALSE) else cat("none\n")
fb <- unique(out$null_used[out$type == "categorical"])
cat("\ncategorical null used:", paste(fb, collapse = ", "), "\n")
