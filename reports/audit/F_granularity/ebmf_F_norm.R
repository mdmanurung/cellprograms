# Amendment A4 real-data sensitivity (pilot-grade, no gate): PFlog / CLR
# normalizations for the singlecell and metacell ebmf arms (workstream F recipe).
# Transforms per 04d-amendment-A4: pflog = zscore_gene(rowcenter(log1p(4*ahat*x/s_c))),
# clr = zscore_gene(rowcenter(log1p(x)))); ahat = NB moment matching on
# depth-normalized counts, median over expressed genes, floor 0.25.
# Concordance vs the completed shiftedlog control fit (F_ebmf_<arm>_loadings.csv).
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellprograms))
suppressPackageStartupMessages(library(Matrix))

root <- "/mnt/shared-workspace/shared/audit/F_granularity"
arms <- c("singlecell", "metacell")
norms <- c("pflog", "clr")
seed <- 20261017

.greedy_match_cor <- function(W1, W2) {
  C <- abs(cor(W1, W2)); C[is.na(C)] <- 0
  total <- 0; n <- min(ncol(W1), ncol(W2))
  for (i in seq_len(n)) {
    idx <- arrayInd(which.max(C), dim(C))
    total <- total + C[idx[1], idx[2]]
    C <- C[-idx[1], -idx[2], drop = FALSE]
    if (nrow(C) == 0 || ncol(C) == 0) break
  }
  total / n
}

.ahat <- function(X) {
  sf <- rowSums(X); sf <- sf / mean(sf)
  Xn <- X / sf
  mu <- colMeans(Xn); v <- apply(Xn, 2, var)
  keep <- mu > 0.01 * max(mu)
  a <- (v[keep] - mu[keep]) / (mu[keep]^2)
  a <- a[is.finite(a) & a > 0]
  max(stats::median(a), 0.25)
}

.transform <- function(X, norm) {
  X <- as.matrix(X)
  if (norm == "pflog") {
    a <- .ahat(X)
    sf <- rowSums(X); sf <- sf / mean(sf)
    Y <- log1p(4 * a * X / sf)
    cat(sprintf("    ahat=%.4f\n", a))
  } else { # clr
    Y <- log1p(X)
  }
  Y <- Y - rowMeans(Y)              # CLR centering (per cell)
  Y <- scale(Y, center = TRUE, scale = TRUE)  # registered gene stage
  Y[, apply(Y, 2, function(j) diff(range(j)) > 0)]
}

rows <- list()
for (arm in arms) {
  d <- file.path(root, "gran", arm)
  X <- as.matrix(readMM(file.path(d, "counts.mtx")))
  meta <- read.delim(file.path(d, "obs_meta.tsv"), stringsAsFactors = FALSE)
  genes <- readLines(file.path(d, "genes.tsv"))
  rownames(X) <- meta$id; colnames(X) <- genes
  ctrl <- read.csv(file.path(root, sprintf("F_ebmf_%s_loadings.csv", arm)), row.names = 1)
  for (norm in norms) {
    cat(sprintf("[%s/%s] obs=%d genes=%d\n", arm, norm, nrow(X), ncol(X)))
    Y <- .transform(X, norm)
    x <- as_cell_program_data(list(A = Y),
         sample_metadata = data.frame(observation_id = meta$id, stringsAsFactors = FALSE))
    t0 <- Sys.time()
    fit <- canonicalize_programs(suppressMessages(suppressWarnings(
      fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = seed))))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    W <- fit$loadings$A
    conc <- .greedy_match_cor(W, ctrl)
    rows[[length(rows) + 1]] <- data.frame(
      method = "ebmf", arm = arm, norm = norm, seed = seed,
      K_selected = ncol(W), runtime_s = round(el, 1),
      concordance_vs_shiftedlog = round(conc, 4))
    write.csv(W, file.path(root, sprintf("F_ebmf_%s_%s_loadings.csv", arm, norm)))
    cat(sprintf("  K=%d %.0fs conc=%.3f\n", ncol(W), el, conc))
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(root, "F_ebmf_norm_sensitivity.csv"), row.names = FALSE)
cat("NORM_SENSITIVITY_DONE\n"); print(out)
