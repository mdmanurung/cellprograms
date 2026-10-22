## run_rep.R -- run all MR-4 arms on one replicate (paired). CLI:
##   Rscript run_rep.R <rep> [outdir]
## Writes <outdir>/rep<rep>.csv (one row per arm) + <outdir>/rep<rep>.rds
## (per-arm fit summaries for the bootstrap in summarize.R).

`%||%` <- function(a, b) if (is.null(a)) b else a
rep <- if (exists("rep") && !is.na(rep)) rep else as.integer(commandArgs(TRUE)[1] %||% "1")
outdir <- commandArgs(TRUE)[2] %||% "out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

source("harness.R")

## registered settings (PREREG_mr4-anchored-v1.md); MR4_G / MR4_DONORS /
## MR4_CELLS / MR4_NMC env overrides exist for the pilot-scale smoke only.
SET <- list(n_donors = as.integer(Sys.getenv("MR4_DONORS", "60")),
            G = as.integer(Sys.getenv("MR4_G", "5000")),
            program_size = 50L, amp = as.numeric(Sys.getenv("MR4_AMP", "1")),
            n_cells_med = as.integer(Sys.getenv("MR4_CELLS", "200")),
            anchor_K = 30L, n_var = as.integer(Sys.getenv("MR4_NVAR", "2000")),
            seed = 4400000L + rep,
            n_perm = as.integer(Sys.getenv("MR4_NPERM", "199")),
            min_donors = as.integer(Sys.getenv("MR4_MINDON", "20")),
            min_shared = as.integer(Sys.getenv("MR4_MINSHARED", "20")),
            n_mc = as.integer(Sys.getenv("MR4_NMC", "20")))

cat(sprintf("[rep %d] simulating (seed %d, amp %.2f)\n", rep, SET$seed, SET$amp))
sim <- sim_cell_level(n_donors = SET$n_donors, G = SET$G, program_size = SET$program_size,
                      amp = SET$amp, n_cells_med = SET$n_cells_med, seed = SET$seed,
                      n_ct = as.integer(Sys.getenv("MR4_NCT", "4")))
X_all <- query_matrix(sim)
cts <- sim$cts
sm_cell <- unique(sim$meta[, c("donor", "site")]); rownames(sm_cell) <- sm_cell$donor
donors <- sm_cell$donor

cat(sprintf("[rep %d] building Stage A anchor (K=%d)\n", rep, SET$anchor_K))
t0 <- proc.time()[["elapsed"]]
anchor <- build_anchor(sim, max_factors = SET$anchor_K, n_var = SET$n_var, seed = 1)
t_anchor <- proc.time()[["elapsed"]] - t0

## donor-mean activities -> caller scores (uniform across arms; means, not
## sums, so donor size does not reweight the correlations)
donor_mean <- function(activities, donor_vec) {
  agg <- rowsum(activities, group = donor_vec)
  agg <- agg[order(rownames(agg)), , drop = FALSE]
  agg / as.numeric(table(donor_vec)[rownames(agg)])
}

run_arms <- list()
## ---- de-novo cell-level comparator ------------------------------------------
cat("[rep] arm denovo_cell\n")
t0 <- proc.time()[["elapsed"]]
sm_cell <- sim$meta[, c("donor", "site")]
sm_cell$observation_id <- rownames(sm_cell)
x_cell <- as_cell_program_data(lapply(cts, function(ct) X_all[sim$meta$cell_type == ct, , drop = FALSE]) |>
                                 stats::setNames(cts),
                               sample_metadata = sm_cell)
fit_dv <- canonicalize_programs(fit_celltype_programs(x_cell, max_factors = 10,
                                                      loading_prior = "point_laplace", var_type = 1,
                                                      features = "variable", n_variable_genes = SET$n_var,
                                                      covariates = "site", covariate_mode = "residualize",
                                                      seed = 1))
t_dv <- proc.time()[["elapsed"]] - t0
run_arms$denovo_cell <- list(fit = fit_dv, t = t_dv)

## ---- anchored arms -----------------------------------------------------------
for (ct in cts) {
  sel <- sim$meta$cell_type == ct
  X_ct <- X_all[sel, , drop = FALSE]
  L <- anchor$loadings[[ct]]                       # genes x K on the anchor panel
  panel <- intersect(rownames(L), colnames(X_ct))
  Xp <- X_ct[, panel, drop = FALSE]
  Lp <- L[panel, , drop = FALSE]
  A0 <- ols_activity(Xp, Lp)
  arms <- list(anch_free = list(fixed_side = "none", add_factors = 0),
               anch_fixed_act = list(fixed_side = "activity", add_factors = 0),
               anch_fixed_eff = list(fixed_side = "effects", add_factors = 0),
               anch_add5 = list(fixed_side = "none", add_factors = 5))
  for (nm in names(arms)) {
    key <- paste0(nm, "_", ct)
    if (is.null(run_arms[[nm]])) run_arms[[nm]] <- list(per_ct = list())
    cat(sprintf("[rep] arm %s / %s\n", nm, ct))
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(arm_anchored(Xp, Lp, A0, fixed_side = arms[[nm]]$fixed_side,
                               add_factors = arms[[nm]]$add_factors, seed = 1),
                  error = function(e) { cat(sprintf("  FAILED %s: %s\n", nm, conditionMessage(e))); NULL })
    if (is.null(r)) next
    run_arms[[nm]]$per_ct[[ct]] <- list(activities = r$activities, loadings = r$effects,
                                        meta = sim$meta[sel, , drop = FALSE])
    run_arms[[nm]]$t <- c(run_arms[[nm]]$t %||% 0, proc.time()[["elapsed"]] - t0)
  }
  ## metacell arm
  key <- "anch_free_mc"
  if (is.null(run_arms[[key]])) run_arms[[key]] <- list(per_ct = list())
  cat(sprintf("[rep] arm anch_free_mc / %s\n", ct))
  t0 <- proc.time()[["elapsed"]]
  g <- metacell_grout(sim$counts[sel, , drop = FALSE], sim$meta[sel, , drop = FALSE],
                      n_mc_per_donor = SET$n_mc, seed = 1)
  cm <- rowsum(sim$counts[sel, , drop = FALSE][g$cell, , drop = FALSE], group = g$mc_id)
  size <- as.numeric(table(g$mc_id)[rownames(cm)])
  Xm <- log1p(t(t(cm) / colSums(cm)) * 1e4)
  Lm <- L[intersect(rownames(L), colnames(Xm)), , drop = FALSE]
  A0m <- ols_activity(Xm[, rownames(Lm)], Lm)
  r <- tryCatch(arm_anchored(Xm[, rownames(Lm)], Lm, A0m, fixed_side = "none",
                             add_factors = 0, seed = 1, S = sqrt(1 / size)),
                error = function(e) { cat(sprintf("  FAILED mc: %s\n", conditionMessage(e))); NULL })
  if (!is.null(r)) {
    mc_donor <- gsub("_mc[0-9]+$", "", rownames(r$activities))
    run_arms[[key]]$per_ct[[ct]] <- list(activities = r$activities, loadings = r$effects,
                                         meta = data.frame(donor = mc_donor, row.names = rownames(r$activities)))
    run_arms[[key]]$t <- c(run_arms[[key]]$t %||% 0, proc.time()[["elapsed"]] - t0)
  }
}

## ---- score every arm: CMP F1 + drift + stability + reconstruction ------------
anchor_scores <- lapply(cts, function(ct) as.matrix(anchor$scores[[ct]])) |> stats::setNames(cts)
truth_pair <- c("A", "B")
rows <- list()
for (nm in names(run_arms)) {
  a <- run_arms[[nm]]
  if (nm == "denovo_cell") {
    scores <- lapply(cts, function(ct) {
      act <- a$fit$scores[[ct]]
      donor_mean(act, donor_of(sim$meta)[rownames(act)])
    }) |> stats::setNames(cts)
    loadings <- a$fit$loadings[cts]
  } else {
    scores <- lapply(cts, function(ct) {
      p <- a$per_ct[[ct]]
      if (is.null(p)) return(NULL)
      donor_mean(p$activities, donor_of(p$meta)[rownames(p$activities)])
    }) |> stats::setNames(cts)
    loadings <- lapply(cts, function(ct) a$per_ct[[ct]]$loadings) |> stats::setNames(cts)
  }
  scores <- scores[!vapply(scores, is.null, logical(1))]
  if (length(scores) < 2L) next
  w <- wrap_for_caller(scores, loadings, donors, sm_cell, names(scores))
  cmp <- tryCatch(find_coordinated_programs(w, ambient = "none", n_perm = SET$n_perm,
                                            min_donors = SET$min_donors, min_shared = SET$min_shared,
                                            seed = 1),
                  error = function(e) { cat(sprintf("  caller FAILED %s: %s\n", nm, conditionMessage(e))); NULL })
  called <- if (is.null(cmp)) NA else paste(cmp$cmps$cell_type, cmp$cmps$program_id, cmp$cmps$cmp,
                                            sep = "/", collapse = ";")
  ## truth: the planted A-B pair is one CMP (same `cmp` id, both types present)
  tp <- if (is.null(cmp)) 0L else {
    ab <- cmp$cmps[cmp$cmps$cell_type %in% truth_pair, ]
    as.integer(length(unique(ab$cmp)) == 1L && nrow(ab) == 2L)
  }
  n_calls <- if (is.null(cmp)) NA_integer_ else length(unique(cmp$cmps$cmp))
  fp <- if (is.null(cmp)) NA_integer_ else n_calls - tp
  f1 <- (2 * tp) / (2 * tp + fp + (1 - tp))
  dr <- mean(vapply(names(scores), function(ct)
    drift_vs_anchor(anchor_scores[[ct]], scores[[ct]]), numeric(1)), na.rm = TRUE)
  st <- mean(vapply(names(scores), function(ct) {
    p <- if (nm == "denovo_cell") NULL else a$per_ct[[ct]]
    L2 <- if (nm == "denovo_cell") a$fit$loadings[[ct]] else p$loadings
    A2 <- if (nm == "denovo_cell") a$fit$scores[[ct]] else p$activities
    ## subset to the fit's own gene panel: ols_activity requires nrow(L) == ncol(X)
    X2 <- if (nm == "denovo_cell") X_all[sim$meta$cell_type == ct, rownames(L2), drop = FALSE] else NULL
    if (nm == "denovo_cell") activity_stability(X2, A2, L2, n_sub = 5L, seed = 1) else NA_real_
  }, numeric(1)), na.rm = TRUE)
  rr <- vapply(names(scores), function(ct) {
    if (nm == "denovo_cell") return(1)
    p <- a$per_ct[[ct]]
    if (is.null(p)) return(NA_real_)
    sel <- sim$meta$cell_type == ct
    Xp <- X_all[sel, rownames(p$loadings), drop = FALSE]
    ref <- run_arms$denovo_cell$fit
    recon_ratio(Xp, p$activities, p$loadings,
                ref$scores[[ct]][rownames(Xp), , drop = FALSE], ref$loadings[[ct]][rownames(p$loadings), , drop = FALSE])
  }, numeric(1))
  rows[[length(rows) + 1L]] <- data.frame(rep = rep, arm = nm, f1 = f1, tp = tp, fp = fp,
                                          n_calls = n_calls, drift = dr, stability = st,
                                          recon_ratio_min = min(rr, na.rm = TRUE),
                                          t_anchor = if (nm == "denovo_cell") t_anchor else NA_real_,
                                          t_arm = if (nm == "denovo_cell") t_dv else sum(a$t %||% 0),
                                          called = called, stringsAsFactors = FALSE)
  cat(sprintf("[rep] %s: F1=%.3f calls=%s drift=%.3f\n", nm, f1, n_calls, dr))
}

out <- do.call(rbind, rows)
write.csv(out, file.path(outdir, sprintf("rep%d.csv", rep)), row.names = FALSE)
saveRDS(list(set = SET, rows = out), file.path(outdir, sprintf("rep%d.rds", rep)))
cat(sprintf("[rep %d] done\n", rep))
