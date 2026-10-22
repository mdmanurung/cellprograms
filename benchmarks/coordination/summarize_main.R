## summarize_main.R — combine per-replicate CSVs into the prereg decision table.
## Paired comparisons (shared seeds): per-replicate recall differences vs pa_cc,
## t-based 95% CI; rates with Wilson 95% CI. run after main_run.R shards finish.
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/summarize_main.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
library(stats)

read_dir <- function(d) {
  fs <- list.files(d, pattern = "\\.csv$", full.names = TRUE)
  if (!length(fs)) return(data.frame())
  do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
}
sub <- Sys.getenv("OUT_SUBDIR", "main")
prog <- read_dir(file.path(here, "out", sub, "prog"))
ds   <- read_dir(file.path(here, "out", sub, "ds"))
cat(sprintf("loaded: %d program rows, %d dataset rows\n", nrow(prog), nrow(ds)))

wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  p <- x / n; z <- 1.96
  den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den)
}

## ---- dataset-level rates per cell family x method ----
ds$any_spurious <- ds$n_spurious > 0
ds$any_trap <- ds$n_trap > 0
fam <- with(ds, paste(mode, breadth, amp, act_cor, sep = "|"))
rate_tab <- do.call(rbind, lapply(split(ds, list(fam, ds$method), drop = TRUE), function(g) {
  n <- nrow(g)
  data.frame(family = paste(g$mode[1], g$breadth[1], g$amp[1], g$act_cor[1], sep = "|"),
             method = g$method[1], n = n,
             P_any_spurious = mean(g$any_spurious),
             mean_spurious = mean(g$n_spurious),
             spurious_per_call = mean(g$spurious_per_call),
             P_any_decoy = mean(g$n_decoy > 0),
             P_any_trap = mean(g$any_trap),
             median_runtime = median(g$runtime_secs),
             stringsAsFactors = FALSE)
}))

## ---- program-level detection per cell family x method ----
fam_p <- with(prog, paste(mode, breadth, amp, act_cor, sep = "|"))
prog_tab <- do.call(rbind, lapply(split(prog, list(fam_p, prog$method), drop = TRUE), function(g) {
  data.frame(family = g$cell_id[1], method = g$method[1],
             program = g$program_id[1],
             n = nrow(g),
             found_rate = mean(g$found),
             mean_recall = mean(g$recall),
             mean_precision = mean(g$precision),
             mean_breadth_error = mean(g$breadth_error, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))

## ---- primary endpoint: paired recall difference vs pa_cc ----
## Primary family (frozen): act_only, each-alone R50 + R20, act_cor = 1 only
## (the 0.9/0.7 levels are the robustness sweep, scored per family in
## summary_programs, never pooled into the primary contrast).
prim <- prog[prog$mode == "act_only" & prog$act_cor == 1 &
               prog$breadth %in% c("R50", "R20"), ]
key <- with(prim, paste(breadth, amp, rep, program_id, sep = "|"))
wide <- lapply(split(prim, prim$method), function(g) setNames(g$recall, with(g, paste(breadth, amp, rep, program_id, sep = "|"))))
stopifnot(all(vapply(wide, function(g) !any(duplicated(names(g))), TRUE)))   # key unique within method
methods <- setdiff(unique(prim$method), "pa_cc")
paired <- if (length(methods) && nrow(prim)) do.call(rbind, lapply(methods, function(m) {
  common <- intersect(names(wide$pa_cc), names(wide[[m]]))
  if (!length(common)) return(NULL)
  dif <- wide[[m]][common] - wide$pa_cc[common]
  data.frame(method = m, n_pairs = length(common),
             mean_recall = mean(wide[[m]][common]),
             pa_cc_recall = mean(wide$pa_cc[common]),
             mean_diff = mean(dif),
             ci_lo = mean(dif) - 1.96 * sd(dif) / sqrt(length(dif)),
             ci_hi = mean(dif) + 1.96 * sd(dif) / sqrt(length(dif)),
             stringsAsFactors = FALSE)
})) else data.frame()
if (is.null(paired)) paired <- data.frame()

write.csv(rate_tab, file.path(here, "out", sub, "summary_rates.csv"), row.names = FALSE)
write.csv(prog_tab, file.path(here, "out", sub, "summary_programs.csv"), row.names = FALSE)
write.csv(paired, file.path(here, "out", sub, "summary_paired_primary.csv"), row.names = FALSE)

cat("\n== primary endpoint: paired recall vs pa_cc (act_only, R50+R20) ==\n")
print(paired, row.names = FALSE)
cat("\n== null family: spurious rates ==\n")
print(rate_tab[grepl("^null", rate_tab$family), ], row.names = FALSE)
cat("\n== priv_corr family: trap rates ==\n")
print(rate_tab[grepl("^priv_corr", rate_tab$family), ], row.names = FALSE)
