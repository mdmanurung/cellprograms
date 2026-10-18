# CP-4 confirm analysis: program_sharing() nuisance/selection sensitivity.
# Evaluates confirm-arm results against PREREGISTRATION.md decision rules
# (incl. Amendment 1). Usage: Rscript analyze_confirm.R <confirm_dir>

args <- commandArgs(trailingOnly = TRUE)
confirm_dir <- args[1]
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
harness_dir <- if (!is.na(file_arg)) dirname(normalizePath(file_arg)) else getwd()

res <- utils::read.csv(file.path(confirm_dir, "results_confirm.csv"), stringsAsFactors = FALSE)
ok <- res[res$status == "ok", , drop = FALSE]
n_err <- sum(res$status != "ok")

ci95 <- function(x, n) {
  if (n == 0) return(c(NA, NA, NA))
  bt <- stats::binom.test(x, n, conf.level = 0.95)
  c(x / n, unname(bt$conf.int[1]), unname(bt$conf.int[2]))
}
f3 <- function(v) if (is.na(v[1])) "—" else sprintf("%.3f [%.3f, %.3f]", v[1], v[2], v[3])

# per-rep promotion indicators
ok$promoted_shared <- ok$classification == "shared"
ok$promoted_any <- ok$n_groups > 0

scen_arm <- unique(ok[, c("scenario", "arm")])
rows <- list()
for (i in seq_len(nrow(scen_arm))) {
  sc <- scen_arm$scenario[i]; ar <- scen_arm$arm[i]
  d <- ok[ok$scenario == sc & ok$arm == ar, , drop = FALSE]
  reps <- unique(d$seed)
  n <- length(reps)
  # per-rep any-shared promotion
  per_rep <- tapply(d$promoted_shared, d$seed, any)
  prom <- ci95(sum(per_rep), n)
  row <- data.frame(scenario = sc, arm = ar, n_reps = n,
                    promo_rate = prom[1], promo_lo = prom[2], promo_hi = prom[3],
                    mean_artifact_frac_top = mean(d$artifact_frac_top[d$promoted_shared]),
                    mean_frac_top_shared = mean(d$frac_top_shared[d$promoted_shared]),
                    gate = NA_character_,
                    stringsAsFactors = FALSE)
  # P1 recovery: 4-type shared group with p_all_shared > 0.5
  if (sc == "P1_shared_biology") {
    rec <- tapply(d$promoted_shared & d$group_size == 4 & d$p_all_shared > 0.5, d$seed, any)
    r <- ci95(sum(rec), n)
    row$gate <- sprintf("recovery %s (gate >= 0.80): %s", f3(r),
                        if (!is.na(r[1]) && r[1] >= 0.8) "PASS" else "INVESTIGATE")
  }
  # N0: promotion rate; defect only if CI excludes 5% high side
  if (sc == "N0_null") {
    row$gate <- sprintf("promotion %s (defect if CI lo > 0.05): %s", f3(prom),
                        if (!is.na(prom[2]) && prom[2] > 0.05) "INVESTIGATE" else "PASS")
  }
  # N1/N3/N4: full-arm promoted groups must carry artifact flag; prespecified must not promote
  if (sc %in% c("N1_vdj_contamination", "N3_small_support", "N4_small_support_vdj")) {
    if (ar == "full") {
      flagged <- mean(d$artifact_frac_top[d$promoted_shared] > 0)
      row$gate <- sprintf("artifact flag on promoted groups = %.2f (gate = 1.00): %s",
                          flagged, if (isTRUE(flagged == 1)) "PASS" else "FLAG GAP")
    } else if (ar == "prespecified") {
      row$gate <- sprintf("prespecified promotion %s (gate = no promotion): %s", f3(prom),
                          if (!is.na(prom[1]) && prom[1] == 0) "PASS" else "INVESTIGATE")
    } else {
      row$gate <- sprintf("posthoc descriptive: promotion %s, mean posthoc_n_drop=%.1f",
                          f3(prom), mean(d$posthoc_n_drop, na.rm = TRUE))
    }
  }
  # N2: full arm must not promote cross-type lineage program
  if (sc == "N2_lineage_vdj") {
    if (ar == "full") {
      row$gate <- sprintf("full-arm promotion %s (gate = 0): %s", f3(prom),
                          if (!is.na(prom[1]) && prom[1] == 0) "PASS" else "INVESTIGATE")
    } else {
      row$gate <- sprintf("%s arm promotion %s (documented cost)", ar, f3(prom))
    }
  }
  rows[[i]] <- row
}
summ <- do.call(rbind, rows)

# Q1 sub-fixture: sharing fires in both regimes; CMP only with correlated scores
q1 <- utils::read.csv(file.path(confirm_dir, "results_q1.csv"), stringsAsFactors = FALSE)
q1_pass <- nrow(q1) == 2 &&
  all(q1$sharing_classification == "shared") &&
  q1$cmp_n[q1$regime == "independent_scores"] == 0 &&
  q1$cmp_n[q1$regime == "correlated_scores"] >= 1

write.csv(summ, file.path(harness_dir, "confirm_summary.csv"), row.names = FALSE)

lines <- c("# CP-4 program-sharing sensitivity — confirmation arm", "",
  sprintf("Date: %s | rows: %d (%d errors) | reps per scenario x arm: %d",
          Sys.Date(), nrow(res), n_err, max(summ$n_reps)), "",
  "| Scenario | Arm | n | shared-promotion [95% CI] | artifact frac (promoted) | Gate |",
  "|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) {
  lines <- c(lines, sprintf("| %s | %s | %d | %s | %.2f | %s |",
    summ$scenario[i], summ$arm[i], summ$n_reps[i],
    f3(c(summ$promo_rate[i], summ$promo_lo[i], summ$promo_hi[i])),
    summ$mean_artifact_frac_top[i], summ$gate[i]))
}
lines <- c(lines, "", "## Q1 two-axis sub-fixture", "",
  sprintf("- sharing fires in both regimes: %s", all(q1$sharing_classification == "shared")),
  sprintf("- CMP fires only with correlated donor scores: %s",
          q1$cmp_n[q1$regime == "independent_scores"] == 0 && q1$cmp_n[q1$regime == "correlated_scores"] >= 1),
  sprintf("- Q1 gate: %s", if (q1_pass) "PASS" else "INVESTIGATE"))
writeLines(lines, file.path(harness_dir, "CONFIRM_ANALYSIS.md"))
cat(paste(lines, collapse = "\n"), "\n")
