## run_b2_b3.R -- B2: V-estimation z_thresh sweep; B3: matching calibration.
## Both run on the UNMASKED K30 fits (existing sharing results are the baseline).
## Outputs: /workspace/v_sensitivity.csv, /workspace/matching_calibration.csv,
## /workspace/b2_b3.rds; checkpointed to /mnt/shared-workspace/cmp/real/.
.libPaths("/workspace/.Rlib")
pkgload::load_all("/workspace/repo2", quiet = TRUE)
stopifnot(requireNamespace("mashr", quietly = TRUE))
R <- "/mnt/shared-workspace/cmp/real"

## ---------- B2: z_thresh sweep ----------
zt_grid <- c(1.5, 2, 2.5, 3)
v_rows <- list()
for (co in c("stephenson", "combat")) {
  fit <- readRDS(file.path(R, sprintf("%s_pb_sex_K30_fit.rds", co)))
  base <- readRDS(sprintf("/workspace/%s_K30_sharing.rds", co))   # z_thresh = 1.5
  base_class <- vapply(base$results, `[[`, "", "classification")
  for (zt in zt_grid) {
    cat("== B2:", co, "z_thresh =", zt, "\n"); flush.console()
    res <- program_sharing(fit, n_perm = 199, seed = 1, z_thresh = zt)
    for (gi in seq_along(res$results)) {
      r <- res$results[[gi]]
      V <- r$V
      off <- V[upper.tri(V)]
      v_rows[[length(v_rows) + 1L]] <- data.frame(
        cohort = co, z_thresh = zt, group_id = gi,
        n_conditions = length(r$cell_types),
        n_null = attr(V, "n_null"), lambda = attr(V, "lambda"),
        v_offdiag_median = if (length(off)) stats::median(off) else NA_real_,
        v_offdiag_max = if (length(off)) max(abs(off)) else NA_real_,
        classification = r$classification,
        flipped_vs_1.5 = if (zt == 1.5) FALSE else r$classification != base_class[gi],
        stringsAsFactors = FALSE)
    }
  }
}
vdf <- do.call(rbind, v_rows)
write.csv(vdf, "/workspace/v_sensitivity.csv", row.names = FALSE)
cat("B2 done. Flip rates by z_thresh:\n")
print(aggregate(flipped_vs_1.5 ~ cohort + z_thresh, vdf, mean))

## ---------- B3: matching calibration at n_perm = 1999 ----------
ari <- function(a, b) {
  # adjusted Rand index for two label vectors (no external deps)
  tab <- table(a, b); n <- sum(tab)
  sum_c <- sum(choose(tab, 2)); sr <- choose(rowSums(tab), 2); sc <- choose(colSums(tab), 2)
  e <- sum(sr) * sum(sc) / choose(n, 2)
  (sum_c - e) / (0.5 * (sum(sr) + sum(sc)) - e)
}
cal_rows <- list()
for (co in c("stephenson", "combat")) {
  cat("== B3:", co, "\n"); flush.console()
  fit <- readRDS(file.path(R, sprintf("%s_pb_sex_K30_fit.rds", co)))
  base <- readRDS(sprintf("/workspace/%s_K30_sharing.rds", co))
  mt199 <- base$matches
  g199 <- base$groups
  t0 <- Sys.time()
  mt2k <- .match_programs_loadings(fit, n_perm = 1999L, seed = 1L)
  cat("matching runtime (n_perm=1999):", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
  g2k <- .sharing_groups(mt2k)

  sig199 <- mt199[mt199$q < 0.05, ]; sig2k <- mt2k[mt2k$q < 0.05, ]
  key <- function(d) paste(d$ct_a, d$prog_a, d$ct_b, d$prog_b)
  e199 <- key(sig199); e2k <- key(sig2k)
  jac <- length(intersect(e199, e2k)) / length(union(e199, e2k))

  # program -> group-id clusterings over the union of matched programs
  memb <- function(groups) {
    m <- list()
    for (i in seq_along(groups)) for (j in seq_len(nrow(groups[[i]])))
      m[[paste(groups[[i]]$cell_type[j], groups[[i]]$program_id[j])]] <- i
    m
  }
  m199 <- memb(g199); m2k <- memb(g2k)
  progs <- union(names(m199), names(m2k))
  lab199 <- vapply(progs, function(p) if (is.null(m199[[p]])) 0L else m199[[p]], 1L)
  lab2k  <- vapply(progs, function(p) if (is.null(m2k[[p]])) 0L else m2k[[p]], 1L)

  cal_rows[[co]] <- data.frame(
    cohort = co,
    n_sig_edges_199 = nrow(sig199), n_sig_edges_1999 = nrow(sig2k),
    edge_jaccard = jac,
    n_groups_199 = length(g199), n_groups_1999 = length(g2k),
    ari_groups = ari(lab199, lab2k),
    cor_q05 = unname(quantile(abs(sig2k$cor), 0.05)),
    cor_q25 = unname(quantile(abs(sig2k$cor), 0.25)),
    cor_q50 = unname(quantile(abs(sig2k$cor), 0.50)),
    cor_q75 = unname(quantile(abs(sig2k$cor), 0.75)),
    stringsAsFactors = FALSE)
}
cal <- do.call(rbind, cal_rows)
write.csv(cal, "/workspace/matching_calibration.csv", row.names = FALSE)
print(cal)
saveRDS(list(v_sensitivity = vdf, matching_calibration = cal), "/workspace/b2_b3.rds")
system2("cp", c("/workspace/b2_b3.rds", file.path(R, "b2_b3.rds")))
cat("ALL DONE\n")
