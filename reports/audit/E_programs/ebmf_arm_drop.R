# Workstream E, arm 1 (dropout detection): cellprograms ebmf pipeline on the
# detection scenario under entrywise count dropout (rates 0.1, 0.3).
# Complements ebmf_arm.R v2 (clean detection, currently running). Same
# shifted-log normalization (vignette recipe, alpha=1, library-size factors).
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellprograms))
out_dir <- "/workspace/audit/cnmf_bench/ebmf"
data_dir <- "/workspace/audit/cnmf_bench/data"

.shifted_log <- function(counts) {
  counts <- as.matrix(counts)
  sf <- rowSums(counts)
  sf <- sf / mean(sf)
  log1p(counts / (1 * sf))
}
.read_counts <- function(f) read.csv(f, row.names = 1, check.names = FALSE)

det <- list()
for (rep in 1:20) for (rate in c(0.1, 0.3)) for (rho in c(0, 0.3, 0.45, 0.6)) {
  d <- readRDS(sprintf("%s/det_drop%g_rho%g_%02d.rds", data_dir, rate, rho, rep))
  slog <- lapply(c("A", "B"), function(ct) {
    .shifted_log(.read_counts(sprintf("%s/det_drop%g_rho%g_%02d_%s_counts.csv", data_dir, rate, rho, rep, ct)))
  })
  names(slog) <- c("A", "B")
  x <- as_cell_program_data(slog, sample_metadata = data.frame(observation_id = rownames(slog$A), site = d$site, stringsAsFactors = FALSE))
  t0 <- Sys.time()
  fit <- suppressMessages(suppressWarnings(fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = rep)))
  fit <- canonicalize_programs(fit)
  res <- suppressWarnings(find_coordinated_programs(fit, strata = "site", n_perm = 199, ambient = "none", seed = rep))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  KA <- ncol(fit$loadings$A); KB <- ncol(fit$loadings$B)
  if (KA == 0 || KB == 0) {
    # at least one type retained no programs (dropout variance inflation) ->
    # cross-type CMPs impossible; skip the permutation call (exact, not approx)
    det[[length(det) + 1]] <- data.frame(arm = "ebmf", rep = rep, rate = rate, rho = rho,
                                         K_A = KA, K_B = KB, detected = FALSE,
                                         false_cmp = FALSE, node_p_A = NA_real_, node_p_B = NA_real_,
                                         runtime_s = el)
    next
  }
  ids_planted <- vapply(c("A", "B"), function(ct) {
    cc <- suppressWarnings(cor(fit$scores[[ct]], d$z))
    colnames(fit$scores[[ct]])[which.max(abs(cc[, 1]))]
  }, "")
  nd <- res$nodes
  pA <- nd$node_p[nd$cell_type == "A" & nd$program_id == ids_planted["A"]]
  pB <- nd$node_p[nd$cell_type == "B" & nd$program_id == ids_planted["B"]]
  cm <- res$cmps
  cA <- cm$cmp[cm$cell_type == "A" & cm$program_id == ids_planted["A"]]
  cB <- cm$cmp[cm$cell_type == "B" & cm$program_id == ids_planted["B"]]
  detected <- length(pA) > 0 && length(pB) > 0 && pA < 0.05 && pB < 0.05 && length(cA) && length(cB) && cA[1] == cB[1]
  det[[length(det) + 1]] <- data.frame(arm = "ebmf", rep = rep, rate = rate, rho = rho,
                                       K_A = KA, K_B = KB, detected = detected,
                                       false_cmp = rho == 0 && nrow(cm) > 0,
                                       node_p_A = pA[1], node_p_B = pB[1], runtime_s = el)
}
det <- do.call(rbind, det)
write.csv(det, sprintf("%s/detection_drop.csv", out_dir), row.names = FALSE)
cat("ebmf dropout detection done. rows:", nrow(det), "\n")
print(aggregate(detected ~ rate + rho, det, FUN = mean))
