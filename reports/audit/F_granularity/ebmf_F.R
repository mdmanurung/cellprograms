# Workstream F (PILOT): ebmf arms on three input granularities.
# Same fit pattern as workstream E v2: shifted-log normalization of COUNTS
# (y = log1p(x/(alpha*s_i)), alpha=1, library-size size factors), then
# fit_celltype_programs(max_factors=10, var_type=1, seed) + canonicalize_programs.
# Three seeds per arm for stability (primary 20261017). Donor scores:
# pseudobulk = scores directly; metacell = n_cells-weighted mean; singlecell = mean.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(cellprograms))
suppressPackageStartupMessages(library(Matrix))

root <- "/mnt/shared-workspace/shared/audit/F_granularity"
arms <- c("pseudobulk", "metacell", "singlecell")
seeds <- c(20261017, 20261018, 20261019)
primary <- as.character(seeds[1])

.shifted_log <- function(counts) {
  counts <- as.matrix(counts)
  sf <- rowSums(counts)
  sf <- sf / mean(sf)
  log1p(counts / (1 * sf))
}

.greedy_match_cor <- function(W1, W2) {
  C <- abs(cor(W1, W2))
  C[is.na(C)] <- 0
  total <- 0
  n <- min(ncol(W1), ncol(W2))
  for (i in seq_len(n)) {
    idx <- arrayInd(which.max(C), dim(C))
    total <- total + C[idx[1], idx[2]]
    C <- C[-idx[1], -idx[2], drop = FALSE]
    if (nrow(C) == 0 || ncol(C) == 0) break
  }
  total / n
}

summary_rows <- list()
for (arm in arms) {
  d <- file.path(root, "gran", arm)
  X <- as.matrix(readMM(file.path(d, "counts.mtx")))
  meta <- read.delim(file.path(d, "obs_meta.tsv"), stringsAsFactors = FALSE)
  genes <- readLines(file.path(d, "genes.tsv"))
  rownames(X) <- meta$id
  colnames(X) <- genes
  cat(sprintf("[%s] obs=%d genes=%d\n", arm, nrow(X), ncol(X)))

  x <- as_cell_program_data(list(A = .shifted_log(X)),
                            sample_metadata = data.frame(observation_id = meta$id, stringsAsFactors = FALSE))
  fits <- list()
  for (s in seeds) {
    t0 <- Sys.time()
    fit <- canonicalize_programs(suppressMessages(suppressWarnings(
      fit_celltype_programs(x, max_factors = 10, var_type = 1, seed = s))))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fits[[as.character(s)]] <- fit
    K <- ncol(fit$loadings$A)
    stab <- if (length(fits) > 1) {
      .greedy_match_cor(fits[[primary]]$loadings$A, fit$loadings$A)
    } else NA_real_
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      method = "ebmf", arm = arm, seed = s, K_selected = K, runtime_s = el,
      seed_stab_loadings = stab)
    cat(sprintf("[%s] seed=%d K=%d %.0fs stab=%s\n", arm, s, K, el,
                ifelse(is.na(stab), "-", sprintf("%.3f", stab))))
  }

  fit <- fits[[primary]]
  W <- fit$loadings$A
  S <- fit$scores$A
  write.csv(W, sprintf("%s/F_ebmf_%s_loadings.csv", root, arm))
  write.csv(S, sprintf("%s/F_ebmf_%s_scores.csv", root, arm))

  # donor-level scores
  if (arm == "pseudobulk") {
    ds <- data.frame(donor = meta$donor, S, check.names = FALSE)
  } else if (arm == "metacell") {
    w <- meta$n_cells
    ds <- do.call(rbind, lapply(split(seq_len(nrow(S)), meta$donor), function(ii) {
      as.data.frame(matrix(colSums(S[ii, , drop = FALSE] * w[ii]),
                           nrow = 1, dimnames = list(NULL, colnames(S))) /
        sum(w[ii]))
    }))
    ds <- data.frame(donor = names(split(seq_len(nrow(S)), meta$donor)), ds, check.names = FALSE)
  } else {
    idx <- split(seq_len(nrow(S)), meta$donor)
    ds <- do.call(rbind, lapply(idx, function(ii) colMeans(S[ii, , drop = FALSE])))
    ds <- data.frame(donor = names(idx), ds, check.names = FALSE)
  }
  write.csv(ds, sprintf("%s/F_ebmf_%s_donor_scores.csv", root, arm), row.names = FALSE)
}
do.call(rbind, summary_rows) |> write.csv(sprintf("%s/F_ebmf_summary.csv", root), row.names = FALSE)
cat("ebmf F arms done\n")
