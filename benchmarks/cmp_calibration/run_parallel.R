# CP-3A calibration — parallel driver for the frozen harness.
#
# Execution detail only: scenario definitions, seeds, n_perm, alpha, tau,
# link_threshold, and decision rules are identical to the serial main() path
# in run_cmp_calibration.R (see PREREGISTRATION.md). Replicates are
# independent and individually seeded, so parallel scheduling does not
# affect results.
#
# Usage:
#   CMP_MC_CORES=12 Rscript run_parallel.R smoke  <out_dir>
#   CMP_MC_CORES=12 Rscript run_parallel.R confirm <out_dir>

suppressPackageStartupMessages(library(cellprograms))

# Source the frozen harness without triggering its main(): harness main()
# only runs when commandArgs(trailingOnly=TRUE) has >= 2 elements at source
# time, so we strip them and restore after.
.args <- commandArgs(trailingOnly = TRUE)
commandArgs <- function(trailingOnly = FALSE) character(0)
source("/workspace/repo/benchmarks/cmp_calibration/run_cmp_calibration.R")
rm(commandArgs)

mode <- match.arg(.args[1], c("smoke", "confirm"))
out_dir <- if (length(.args) >= 2) .args[2] else file.path("/workspace/cmp_calib", mode)
n_reps <- if (mode == "smoke") 20L else 500L
n_perm <- if (mode == "smoke") 199L else 499L
mc_cores <- as.integer(Sys.getenv("CMP_MC_CORES", "12"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scns <- .scenarios()
t_start <- Sys.time()
message(sprintf("[%s] start %s | %d scenarios x %d reps | n_perm=%d | cores=%d",
                mode, format(t_start, tz = "UTC", usetz = TRUE),
                length(scns), n_reps, n_perm, mc_cores))

.to_df <- function(x) {
  if (inherits(x, "try-error")) {
    return(data.frame(scenario = NA_character_, seed = NA_integer_,
                      status = "error", message = as.character(x),
                      elapsed = NA_real_, n_nodes = NA, n_sig = NA, any_sig = NA,
                      null_n_sig = NA, null_n = NA, plant_n_sig = NA, plant_n = NA,
                      n_cmps = NA, cmp_membership = "", cmp_signs = "",
                      stringsAsFactors = FALSE))
  }
  if (identical(x$status, "error")) {
    return(data.frame(scenario = x$scenario, seed = x$seed, status = "error",
                      message = x$message, elapsed = x$elapsed,
                      n_nodes = NA, n_sig = NA, any_sig = NA,
                      null_n_sig = NA, null_n = NA, plant_n_sig = NA, plant_n = NA,
                      n_cmps = NA, cmp_membership = "", cmp_signs = "",
                      stringsAsFactors = FALSE))
  }
  data.frame(scenario = x$scenario, seed = x$seed, status = x$status,
             message = "", elapsed = x$elapsed, n_nodes = x$n_nodes,
             n_sig = x$n_sig, any_sig = x$any_sig,
             null_n_sig = x$null_n_sig, null_n = x$null_n,
             plant_n_sig = x$plant_n_sig, plant_n = x$plant_n,
             n_cmps = x$n_cmps, cmp_membership = x$cmp_membership,
             cmp_signs = x$cmp_signs, stringsAsFactors = FALSE)
}

all_df <- list()
for (si in seq_along(scns)) {
  t0 <- Sys.time()
  res_list <- parallel::mclapply(seq_len(n_reps), function(r) {
    seed <- if (mode == "smoke") 9000L + r else 202610180L + si * 1000L + r
    tryCatch(run_replicate(names(scns)[si], scns[[si]], seed, n_perm),
             error = function(e) structure(conditionMessage(e), class = "try-error"))
  }, mc.cores = mc_cores, mc.preschedule = TRUE)
  df <- do.call(rbind, lapply(res_list, .to_df))
  all_df[[si]] <- df
  saveRDS(res_list, file.path(out_dir, sprintf("checkpoint_%s_%02d_%s.rds",
                                               mode, si, names(scns)[si])))
  write.csv(df, file.path(out_dir, sprintf("results_%s_%02d_%s.csv",
                                           mode, si, names(scns)[si])), row.names = FALSE)
  message(sprintf("[%s] %s done: %d rows, %d errors, %.1f min",
                  mode, names(scns)[si], nrow(df), sum(df$status == "error"),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

df <- do.call(rbind, all_df)
write.csv(df, file.path(out_dir, paste0("results_", mode, ".csv")), row.names = FALSE)

manifest <- list(mode = mode, n_reps = n_reps, n_perm = n_perm, alpha = 0.05,
                 tau = 2, link_threshold = 3,
                 scenarios = names(scns),
                 seed_rule = if (mode == "smoke") "9000 + replicate" else
                   "202610180 + scenario*1000 + replicate",
                 driver = "run_parallel.R",
                 mc_cores = mc_cores,
                 cellprograms_sha = system("git -C /workspace/repo rev-parse HEAD", intern = TRUE),
                 r_version = R.version.string,
                 started = format(t_start, tz = "UTC", usetz = TRUE),
                 finished = format(Sys.time(), tz = "UTC", usetz = TRUE),
                 wall_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "run_manifest.json"))
message(sprintf("[%s] wrote %d replicate rows to %s in %.1f min",
                mode, nrow(df), out_dir, manifest$wall_minutes))
