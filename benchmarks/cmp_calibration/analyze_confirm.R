# CP-3A confirmation-arm analysis: per-scenario operating characteristics with
# exact binomial CIs, evaluated against the preregistered decision rules
# (PREREGISTRATION.md). Inputs: results_confirm_*.csv in confirm_dir.
# Outputs: confirm_summary.csv + CONFIRM_ANALYSIS.md in the harness dir.
#
# Usage: Rscript analyze_confirm.R <confirm_dir> [interrupted_dir]

args <- commandArgs(trailingOnly = TRUE)
confirm_dir <- args[1]
interrupted_dir <- if (length(args) >= 2) args[2] else NULL
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
harness_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else getwd()

files <- list.files(confirm_dir, pattern = "^results_confirm_[0-9]+_.*\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
res <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
res <- res[res$status == "ok", , drop = FALSE]

alpha <- 0.05
binom_ci <- function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_, NA_real_))
  bt <- stats::binom.test(x, n, conf.level = 0.95)
  c(x / n, unname(bt$conf.int[1]), unname(bt$conf.int[2]))
}

scenarios <- sort(unique(res$scenario))
rows <- list()
for (sc in scenarios) {
  d <- res[res$scenario == sc, , drop = FALSE]
  n <- nrow(d)
  fam <- sub("_.*$", "", sc)
  row <- data.frame(scenario = sc, n_valid = n,
                    fwer = NA_real_, fwer_lo = NA_real_, fwer_hi = NA_real_,
                    null_fpr = NA_real_, null_fpr_lo = NA_real_, null_fpr_hi = NA_real_,
                    plant_power = NA_real_, plant_power_lo = NA_real_, plant_power_hi = NA_real_,
                    mean_cmps = NA_real_, sign_mode = NA_character_, sign_mode_frac = NA_real_,
                    stringsAsFactors = FALSE)
  # family-wise: any significant CMP node
  fwer <- binom_ci(sum(d$any_sig), n)
  row$fwer <- fwer[1]; row$fwer_lo <- fwer[2]; row$fwer_hi <- fwer[3]
  # null-node FPR (B1): per-null-node rate
  if (any(!is.na(d$null_n) & d$null_n > 0)) {
    nf <- binom_ci(sum(d$null_n_sig, na.rm = TRUE), sum(d$null_n, na.rm = TRUE))
    row$null_fpr <- nf[1]; row$null_fpr_lo <- nf[2]; row$null_fpr_hi <- nf[3]
  }
  # planted-node recovery (power) where plants exist
  if (any(!is.na(d$plant_n) & d$plant_n > 0)) {
    pw <- binom_ci(sum(d$plant_n_sig, na.rm = TRUE), sum(d$plant_n, na.rm = TRUE))
    row$plant_power <- pw[1]; row$plant_power_lo <- pw[2]; row$plant_power_hi <- pw[3]
  }
  row$mean_cmps <- mean(d$n_cmps)
  if (fam == "D") {
    signs <- table(d$cmp_signs)
    row$sign_mode <- names(sort(signs, decreasing = TRUE))[1]
    row$sign_mode_frac <- max(signs) / n
  }
  rows[[sc]] <- row
}
summ <- do.call(rbind, rows)
rownames(summ) <- NULL

# decision rules
summ$family <- sub("_.*$", "", summ$scenario)
summ$verdict <- NA_character_
for (i in seq_len(nrow(summ))) {
  fam <- summ$family[i]
  if (fam %in% c("A0", "A1")) {
    covers <- summ$fwer_lo[i] <= alpha && summ$fwer_hi[i] >= alpha
    summ$verdict[i] <- if (covers) "CALIBRATED (CI covers alpha)" else "MISCALIBRATED (CI excludes alpha)"
  } else if (fam == "B1") {
    summ$verdict[i] <- if (!is.na(summ$null_fpr_hi[i]) && summ$null_fpr_hi[i] <= alpha)
      "STRONG FWER CONTROL (null-node FPR CI <= alpha)" else "STRONG FWER CLAIM NOT SUPPORTED"
  } else if (fam %in% c("C1", "C2", "C3", "C4")) {
    summ$verdict[i] <- sprintf("FWER estimate %.3f [%.3f, %.3f] vs alpha %.2f -> document",
                               summ$fwer[i], summ$fwer_lo[i], summ$fwer_hi[i], alpha)
  } else if (fam %in% c("D1", "D2")) {
    summ$verdict[i] <- "DESCRIPTIVE"
  }
}

# determinism check vs interrupted run (same frozen seeds, A0-C3)
determinism <- NA_character_
if (!is.null(interrupted_dir) && dir.exists(interrupted_dir)) {
  ifiles <- list.files(interrupted_dir, pattern = "^results_confirm_[0-9]+_.*\\.csv$", full.names = TRUE)
  if (length(ifiles) > 0) {
    old <- do.call(rbind, lapply(ifiles, utils::read.csv, stringsAsFactors = FALSE))
    old <- old[old$status == "ok", c("scenario", "seed", "any_sig", "n_sig")]
    new <- res[, c("scenario", "seed", "any_sig", "n_sig")]
    m <- merge(old, new, by = c("scenario", "seed"), suffixes = c("_old", "_new"))
    n_match <- sum(m$any_sig_old == m$any_sig_new & m$n_sig_old == m$n_sig_new)
    determinism <- sprintf("%d/%d replicates identical (any_sig, n_sig)", n_match, nrow(m))
  }
}

write.csv(summ, file.path(harness_dir, "confirm_summary.csv"), row.names = FALSE)

lines <- c(
  "# CP-3A CMP calibration — confirmation arm",
  "",
  sprintf("Date: %s | replicates/scenario: %d | n_perm = 499 | alpha = 0.05",
          Sys.Date(), max(summ$n_valid)),
  sprintf("Valid rows: %d | Determinism vs interrupted run: %s", nrow(res), determinism),
  "",
  "| Scenario | n | FWER [95% CI] | null-node FPR [95% CI] | plant power [95% CI] | Verdict |",
  "|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(summ))) {
  fmt <- function(est, lo, hi) if (is.na(est)) "—" else sprintf("%.3f [%.3f, %.3f]", est, lo, hi)
  lines <- c(lines, sprintf("| %s | %d | %s | %s | %s | %s |",
    summ$scenario[i], summ$n_valid[i],
    fmt(summ$fwer[i], summ$fwer_lo[i], summ$fwer_hi[i]),
    fmt(summ$null_fpr[i], summ$null_fpr_lo[i], summ$null_fpr_hi[i]),
    fmt(summ$plant_power[i], summ$plant_power_lo[i], summ$plant_power_hi[i]),
    summ$verdict[i]))
}
if (any(grepl("^D", summ$family))) {
  lines <- c(lines, "", "## D-family sign patterns", "")
  for (i in which(summ$family == "D")) {
    lines <- c(lines, sprintf("- %s: modal cmp_signs = %s (%.1f%% of replicates)",
                              summ$scenario[i], summ$sign_mode[i], 100 * summ$sign_mode_frac[i]))
  }
}
writeLines(lines, file.path(harness_dir, "CONFIRM_ANALYSIS.md"))
cat(paste(lines, collapse = "\n"), "\n")
