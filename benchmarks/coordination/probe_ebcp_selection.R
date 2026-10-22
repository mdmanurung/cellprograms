## probe_ebcp_selection.R — P3b diagnosis: why does the eb_cp arm lose power on
## act+genes?  For each act+genes R50 replicate (seeds from the frozen main
## grid), sweep the two selection knobs (Kmax, heldout_z) and record, per
## greedy component: held-out z (trace), whether kept, and the correlation of
## the refit donor activity (V) with the true shared z and with the nuisance
## institute z.  Answers: (1) does the true program's component clear the
## frozen heldout_z=3 gate; (2) is power lost at the z gate or at membership;
## (3) does Kmax bind before the z gate stops the greedy pass.
## run: R_LIBS_USER=/workspace/.Rlib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
##      Rscript benchmarks/coordination/probe_ebcp_selection.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "score.R")) source(file.path(here, f))

seeds <- c(2222L, 2260L, 2298L)          # act+genes R50 a0.5 c1 reps 1-3 (main grid)
Kmaxs <- c(15L, 30L)
hzs   <- c(3, 2, 1)
options(ms.G = 2000L)

rows <- list()
for (rep_i in seq_along(seeds)) {
  seed <- seeds[rep_i]
  d <- gen_coord("act+genes", amp = 0.5, seed = seed, act_cor = 1, breadth = "R50")
  fit <- canonicalize_programs(fit_celltype_programs(
    d$x, max_factors = 10L, covariates = "institute", covariate_mode = "residualize", seed = seed))
  mats <- lapply(fit$fits[fit$cell_types], function(f) f$Y_used)
  z_true <- d$truth$programs$R50$z[[1]]             # shared z (act_cor = 1), named by donor
  bz     <- d$truth$nuisance
  mem    <- d$truth$programs$R50$members
  for (K in Kmaxs) for (hz in hzs) {
    t0 <- Sys.time()
    cp <- eb_cp_fit(mats, Kmax = K, heldout_z = hz, strata = d$strata, seed = 1L)
    t_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    n_comp <- cp$n_components
    for (r in seq_len(n_comp)) {
      v <- cp$V[, r]
      cor_true <- suppressWarnings(cor(v, z_true[rownames(cp$V)]))
      cor_bz   <- suppressWarnings(cor(v, bz[rownames(cp$V)]))
      ## membership at the frozen threshold: |W| max-normalized within component
      w_rel <- abs(cp$W[, r]) / max(abs(cp$W[, r]))
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_i, seed = seed, Kmax = K, heldout_z = hz,
        component = r, n_components = n_comp,
        kept_z = cp$kept_heldout_z[r], trace_z = cp$trace_z[r],
        cor_true_z = round(cor_true, 3), cor_nuisance = round(cor_bz, 3),
        n_mem = sum(w_rel >= cp$member_thresh), n_true_mem = sum(w_rel[mem] >= cp$member_thresh),
        is_best_true = r == which.max(abs(apply(cp$V, 2, function(vv)
          suppressWarnings(cor(vv, z_true[rownames(cp$V)]))))),
        secs = round(t_s, 1), stringsAsFactors = FALSE)
    }
    if (n_comp == 0L) rows[[length(rows) + 1L]] <- data.frame(
      rep = rep_i, seed = seed, Kmax = K, heldout_z = hz, component = NA_integer_,
      n_components = 0L, kept_z = NA_real_, trace_z = if (length(cp$trace_z)) cp$trace_z[1] else NA_real_,
      cor_true_z = NA_real_, cor_nuisance = NA_real_, n_mem = NA_integer_,
      n_true_mem = NA_integer_, is_best_true = NA, secs = round(t_s, 1),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("[rep %d seed %d] done\n", rep_i, seed))
}
out <- do.call(rbind, rows)
write.csv(out, file.path(here, "out", "probe_ebcp_selection.csv"), row.names = FALSE)

## compact per-config view: first-component trace_z, n_components, best true-z cor
agg <- do.call(rbind, lapply(split(out, paste(out$rep, out$Kmax, out$heldout_z)), function(g) {
  data.frame(rep = g$rep[1], Kmax = g$Kmax[1], heldout_z = g$heldout_z[1],
             n_comp = g$n_components[1],
             trace_z_c1 = round(g$trace_z[1], 2),
             trace_z_c2 = if (nrow(g) > 1) round(g$trace_z[2], 2) else NA_real_,
             best_cor_true = round(max(abs(g$cor_true_z), na.rm = TRUE), 3),
             best_comp = g$component[which.max(abs(g$cor_true_z))],
             n_true_mem_best = g$n_true_mem[which.max(abs(g$cor_true_z))],
             stringsAsFactors = FALSE)
}))
cat("\n=== eb_cp selection diagnosis (act+genes R50, amp 0.5, act_cor 1) ===\n")
print(agg, row.names = FALSE)
cat("\nFull per-component table: out/probe_ebcp_selection.csv\n")
