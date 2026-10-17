# Workstream E, arm 1 v2: cellprograms ebmf pipeline on the shared benchmark datasets.
# v2: input normalization fixed to the flashier single-cell vignette recipe —
# shifted log counts y = log(1 + x/(alpha*s_i)), alpha = 1, library-size size
# factors — applied to the COUNT matrices (same input scale as the cNMF arm).
# Adds dropout scenarios (drop0.1, drop0.3). The vignette's variance floor
# (S = sd(log(Poisson(1/n)+1))) is NOT added: fit_celltype_programs does not
# expose S; that is CAND-13, to be tested separately.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellprograms))
out_dir <- "/workspace/audit/cnmf_bench/ebmf"
dir.create(out_dir, showWarnings = FALSE)
data_dir <- "/workspace/audit/cnmf_bench/data"

.shifted_log <- function(counts) {
  counts <- as.matrix(counts)
  sf <- rowSums(counts)
  sf <- sf / mean(sf)
  log1p(counts / (1 * sf))
}
.read_counts <- function(f) read.csv(f, row.names = 1, check.names = FALSE)
.top_jaccard <- function(W, block_genes, top = 30) {
  vapply(seq_len(ncol(W)), function(j) {
    g <- head(rownames(W)[order(-abs(W[, j]))], top)
    length(intersect(g, block_genes)) / length(union(g, block_genes))
  }, 1)
}

# ---------- recovery ----------
scenarios <- list(
  list(tag = "k3", reps = 1:10, k = 3L),
  list(tag = "noise", reps = 11:20, k = 0L),
  list(tag = "drop0.1", reps = 1:10, k = 3L),
  list(tag = "drop0.3", reps = 1:10, k = 3L)
)
rec <- list()
for (sc in scenarios) {
  for (rep in sc$reps) {
    counts <- .read_counts(sprintf("%s/recovery_%s_%02d_counts.csv", data_dir, sc$tag, rep))
    d <- readRDS(sprintf("%s/recovery_%s_%02d.rds", data_dir, sc$tag, rep))
    x <- as_cell_program_data(list(A = .shifted_log(counts)),
                              sample_metadata = data.frame(observation_id = rownames(counts), stringsAsFactors = FALSE))
    t0 <- Sys.time()
    fit <- suppressMessages(suppressWarnings(
      fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = rep)))
    fit <- canonicalize_programs(fit)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    K <- ncol(fit$loadings$A)
    S <- fit$scores$A
    if (d$k_true > 0 && K > 0) {
      cc <- suppressWarnings(cor(S, d$z))
      score_cor <- apply(cc, 2, max)
      W <- fit$loadings$A
      jac <- vapply(seq_len(d$k_true), function(k) {
        bg <- rownames(d$W)[d$W[, k] == 1]
        max(.top_jaccard(W, bg))
      }, 1)
    } else { score_cor <- NA_real_; jac <- NA_real_ }
    rec[[length(rec) + 1]] <- data.frame(arm = "ebmf", scenario = sc$tag, rep = rep, K_selected = K,
                                         k_true = d$k_true,
                                         score_cor = if (d$k_true > 0) mean(score_cor, na.rm = TRUE) else NA_real_,
                                         top_jaccard = if (d$k_true > 0) mean(jac, na.rm = TRUE) else NA_real_,
                                         runtime_s = el)
  }
}
# seed reproducibility: 2 pairs on recovery k3 rep 1 (shifted-log input)
counts <- .read_counts(sprintf("%s/recovery_k3_01_counts.csv", data_dir))
x <- as_cell_program_data(list(A = .shifted_log(counts)),
                          sample_metadata = data.frame(observation_id = rownames(counts), stringsAsFactors = FALSE))
for (s in c(1, 2)) {
  f1 <- canonicalize_programs(suppressMessages(suppressWarnings(fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = s))))
  f2 <- canonicalize_programs(suppressMessages(suppressWarnings(fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = s + 100))))
  rec[[length(rec) + 1]] <- data.frame(arm = "ebmf", scenario = "seedpair", rep = s,
                                       K_selected = ncol(f1$loadings$A), k_true = 3,
                                       score_cor = max(suppressWarnings(cor(f1$scores$A, f2$scores$A))),
                                       top_jaccard = NA, runtime_s = NA)
}
rec <- do.call(rbind, rec)
write.csv(rec, sprintf("%s/recovery.csv", out_dir), row.names = FALSE)

# ---------- detection ----------
det <- list()
for (rep in 1:20) for (rho in c(0, 0.3, 0.45, 0.6)) {
  d <- readRDS(sprintf("%s/det_rho%g_%02d.rds", data_dir, rho, rep))
  slog <- lapply(c("A", "B"), function(ct) {
    .shifted_log(.read_counts(sprintf("%s/det_rho%g_%02d_%s_counts.csv", data_dir, rho, rep, ct)))
  })
  names(slog) <- c("A", "B")
  x <- as_cell_program_data(slog, sample_metadata = data.frame(observation_id = rownames(slog$A), site = d$site, stringsAsFactors = FALSE))
  t0 <- Sys.time()
  fit <- suppressMessages(suppressWarnings(fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = rep)))
  fit <- canonicalize_programs(fit)
  res <- suppressWarnings(find_coordinated_programs(fit, strata = "site", n_perm = 199, ambient = "none", seed = rep))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
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
  det[[length(det) + 1]] <- data.frame(arm = "ebmf", rep = rep, rho = rho, detected = detected,
                                       false_cmp = rho == 0 && nrow(cm) > 0,
                                       node_p_A = pA[1], node_p_B = pB[1], runtime_s = el)
}
det <- do.call(rbind, det)
write.csv(det, sprintf("%s/detection.csv", out_dir), row.names = FALSE)
cat("ebmf arm v2 done. recovery rows:", nrow(rec), " detection rows:", nrow(det), "\n")
print(aggregate(detected ~ rho, det, FUN = mean))
