## probe_blockfix.R — validate the permutation-calibrated block_E on one dataset.
## env: MODE, BREADTH, REP, SEED, OUT (tag for the log line)
## run: MODE=null BREADTH=R50 REP=1 SEED=... OUT=null R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/probe_blockfix.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
for (f in c("sim_coord.R", "methods.R", "score.R")) source(file.path(here, f))
mode <- Sys.getenv("MODE"); breadth <- Sys.getenv("BREADTH")
rep_i <- as.integer(Sys.getenv("REP", "1")); seed <- as.integer(Sys.getenv("SEED"))
tag <- Sys.getenv("OUT", "probe")
options(ms.G = 2000L)

d <- gen_coord(mode, amp = 0.5, seed = seed, breadth = breadth)
fit <- canonicalize_programs(fit_celltype_programs(
  d$x, max_factors = 10L, covariates = "institute", covariate_mode = "residualize", seed = seed))

t0 <- Sys.time()
raw <- method_block_e(fit)
t_raw <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
s_raw <- score_method(raw, d$truth)

t0 <- Sys.time()
cal <- method_block_e(fit, strata = d$strata, n_perm = 199L, alpha = 0.05)
t_cal <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
s_cal <- score_method(cal, d$truth)

cat(sprintf("[%s] %s r%03d seed %d\n", tag, mode, rep_i, seed))
cat(sprintf("  raw : n_calls=%d n_spurious=%d n_trap=%d found=%s (%.0f s)\n",
            s_raw$dataset_level$n_calls, s_raw$dataset_level$n_spurious,
            s_raw$dataset_level$n_trap,
            paste(s_raw$program_level$found, collapse = ","), t_raw))
cat(sprintf("  cal : n_calls=%d n_spurious=%d n_trap=%d found=%s (%.0f s)\n",
            s_cal$dataset_level$n_calls, s_cal$dataset_level$n_spurious,
            s_cal$dataset_level$n_trap,
            paste(s_cal$program_level$found, collapse = ","), t_cal))
if (!is.null(cal$factor_p)) cat("  factor p:", paste(sprintf("%s=%.3f", cal$factor_p$factor_id, cal$factor_p$p), collapse = " "), "\n")
