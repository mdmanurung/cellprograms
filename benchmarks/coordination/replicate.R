## replicate.R — one main-run replicate: dataset -> stage-1 (once) -> all arms
## -> scorer. Per-replicate CSVs in <outdir>/{prog,ds}/ make runs resumable.
## Library file: caller sources sim_coord.R / methods.R / eb_cp.R / score.R and
## the package R/ files first. Driver: main_run.R.

METHODS_ALL <- c("pa_cc", "node_cor", "subspace_sum", "gene_match", "block_E", "eb_cp",
                 "oracle")

run_one_replicate <- function(mode, breadth, amp, rep, seed, act_cor = 1,
                              G = 2000L, outdir, methods = METHODS_ALL) {
  options(ms.G = as.integer(G))
  cell_id <- sprintf("%s_%s_a%g_c%g", mode, breadth, amp, act_cor)
  fprog <- file.path(outdir, "prog", sprintf("%s_r%03d.csv", cell_id, rep))
  fds   <- file.path(outdir, "ds",   sprintf("%s_r%03d.csv", cell_id, rep))
  if (file.exists(fds)) return(invisible("skip"))
  d <- gen_coord(mode, amp = amp, seed = seed, act_cor = act_cor, breadth = breadth)

  t0 <- Sys.time()
  fit <- fit_celltype_programs(d$x, max_factors = 10L, covariates = "institute",
                               covariate_mode = "residualize", seed = seed)
  fit <- canonicalize_programs(fit)
  t_stage1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  prog_rows <- list(); ds_rows <- list()
  for (m in methods) {
    t0 <- Sys.time()
    res <- tryCatch(
      if (m == "oracle") {
        method_oracle(lapply(fit$fits[fit$cell_types], function(f) f$Y_used),
                      d$truth, d$strata, n_perm = 999L, seed = seed)
      } else {
        run_method(m, fit, d$strata, seed = seed)
      }, error = function(e) e)
    t_m <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
      ds_rows[[m]] <- data.frame(cell_id = cell_id, mode = mode, breadth = breadth,
                                 amp = amp, act_cor = act_cor, rep = rep, seed = seed,
                                 method = m, error = conditionMessage(res),
                                 runtime_secs = t_m, stage1_secs = t_stage1,
                                 n_calls = NA_integer_, n_spurious = NA_integer_,
                                 n_decoy = NA_integer_, n_trap = NA_integer_,
                                 spurious_per_call = NA_real_, stringsAsFactors = FALSE)
      next
    }
    s <- score_method(res, d$truth)
    dl <- cbind(data.frame(cell_id = cell_id, mode = mode, breadth = breadth,
                           amp = amp, act_cor = act_cor, rep = rep, seed = seed,
                           error = NA_character_, runtime_secs = t_m,
                           stage1_secs = t_stage1, stringsAsFactors = FALSE),
                s$dataset_level)
    ds_rows[[m]] <- dl
    if (nrow(s$program_level)) {
      prog_rows[[m]] <- cbind(dl[, c("cell_id", "mode", "breadth", "amp",
                                     "act_cor", "rep", "seed", "method")],
                              s$program_level)
    }
  }
  prog_df <- if (length(prog_rows)) do.call(rbind, prog_rows) else data.frame(
    cell_id = character(), mode = character(), breadth = character(),
    amp = numeric(), act_cor = numeric(), rep = integer(), seed = integer(),
    method = character(), program_id = character(), n_true = integer(),
    matched_cluster = integer(), recall = numeric(), precision = numeric(),
    breadth_error = integer(), found = logical(), stringsAsFactors = FALSE)
  write.csv(prog_df, fprog, row.names = FALSE)
  write.csv(do.call(rbind, ds_rows), fds, row.names = FALSE)
  invisible("done")
}
