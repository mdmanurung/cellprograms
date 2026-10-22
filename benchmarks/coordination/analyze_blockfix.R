## analyze_blockfix.R — summary for the calibrated block_E rerun (out/blockfix)
## vs the frozen main run (out/main, identical seeds). Gates (decision-table
## follow-up): null FPR ~ 0 (raw: 1.00), priv_corr trap ~ 0 (raw: 1.00),
## primary recall gain vs pa_cc preserved (raw block_E: +0.130 [+0.004, +0.255]).
## Label: post-prereg development iteration evaluated against the same frozen
## gates — NOT confirmatory.
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
library(stats)

read_dir <- function(d) {
  fs <- list.files(d, pattern = "\\.csv$", full.names = TRUE)
  if (!length(fs)) return(data.frame())
  do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
}
cal_ds <- read_dir(file.path(here, "out", "blockfix", "ds"))
cal_pg <- read_dir(file.path(here, "out", "blockfix", "prog"))
main_ds <- read_dir(file.path(here, "out", "main", "ds"))
main_pg <- read_dir(file.path(here, "out", "main", "prog"))
cat(sprintf("loaded: blockfix %d ds / %d prog rows; main %d ds / %d prog rows\n",
            nrow(cal_ds), nrow(cal_pg), nrow(main_ds), nrow(main_pg)))
stopifnot(nrow(cal_ds) > 0)

wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  p <- x / n; z <- 1.96; den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den)
}
rate_row <- function(g, label) {
  n <- nrow(g)
  w1 <- wilson(sum(g$n_spurious > 0), n); w2 <- wilson(sum(g$n_trap > 0), n)
  data.frame(arm = label, family = paste(g$mode[1], g$breadth[1], g$amp[1], g$act_cor[1], sep = "|"),
             n = n, P_any_spurious = mean(g$n_spurious > 0), ci_lo_spur = w1[1], ci_hi_spur = w1[2],
             mean_spurious = round(mean(g$n_spurious), 1),
             P_any_trap = mean(g$n_trap > 0), ci_lo_trap = w2[1], ci_hi_trap = w2[2],
             median_runtime = round(median(g$runtime_secs), 1), stringsAsFactors = FALSE)
}
fam_c <- with(cal_ds, paste(mode, breadth, amp, act_cor, sep = "|"))
rates <- do.call(rbind, c(
  lapply(split(cal_ds, fam_c), rate_row, label = "block_E_cal"),
  lapply(split(main_ds[main_ds$method == "block_E", ],
               with(main_ds[main_ds$method == "block_E", ], paste(mode, breadth, amp, act_cor, sep = "|"))),
         rate_row, label = "block_E_raw"),
  lapply(split(main_ds[main_ds$method == "pa_cc", ],
               with(main_ds[main_ds$method == "pa_cc", ], paste(mode, breadth, amp, act_cor, sep = "|"))),
         rate_row, label = "pa_cc")))
rates <- rates[order(rates$family, rates$arm), ]
cat("\n== dataset-level rates (calibrated vs raw vs pa_cc) ==\n")
print(rates, row.names = FALSE)

## paired program-level contrasts on the primary family (act_only R50+R20, act_cor 1)
paired_contrast <- function(pg, m, base = "pa_cc") {
  prim <- pg[pg$mode == "act_only" & pg$act_cor == 1 & pg$breadth %in% c("R50", "R20"), ]
  wide <- lapply(split(prim, prim$method), function(g)
    setNames(g$recall, with(g, paste(breadth, amp, rep, program_id, sep = "|"))))
  if (!all(c(base, m) %in% names(wide))) return(NULL)
  common <- intersect(names(wide[[base]]), names(wide[[m]]))
  if (!length(common)) return(NULL)
  dif <- wide[[m]][common] - wide[[base]][common]
  data.frame(method = m, vs = base, n_pairs = length(common),
             mean_recall = round(mean(wide[[m]][common]), 3),
             base_recall = round(mean(wide[[base]][common]), 3),
             mean_diff = round(mean(dif), 3),
             ci_lo = round(mean(dif) - 1.96 * sd(dif) / sqrt(length(dif)), 3),
             ci_hi = round(mean(dif) + 1.96 * sd(dif) / sqrt(length(dif)), 3),
             stringsAsFactors = FALSE)
}
pc <- do.call(rbind, c(paired_contrast(cal_pg, "block_E_cal", "pa_cc"),
                       paired_contrast(main_pg, "block_E", "pa_cc"),
                       paired_contrast(cal_pg, "block_E_cal", "block_E")))
cat("\n== paired recall contrasts (primary family) ==\n")
print(pc, row.names = FALSE)

write.csv(rates, file.path(here, "out", "blockfix_summary_rates.csv"), row.names = FALSE)
if (!is.null(pc)) write.csv(pc, file.path(here, "out", "blockfix_summary_paired.csv"), row.names = FALSE)
