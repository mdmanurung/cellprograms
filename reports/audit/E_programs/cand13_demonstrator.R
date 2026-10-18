# CAND-13 demonstrator: variance-floor S calibration (registered amendment 3)
#
# Question (ledger 2026-10-17, amendment 3): does the flashier single-cell
# vignette's variance lower bound S = sd(log(Poisson(1/n)+1)) improve lfsr
# calibration, and does it mitigate the entrywise-dropout collapse (which is
# driven by too-LARGE variance estimates)? The vignette passes S as a FIXED
# residual SD (flashier stops estimating variances when S is supplied), so S
# can in principle prevent both miscalibration directions. Both are tested.
#
# Part A: planted-truth dropout scenarios (E-workstream recovery_k3 generator,
#         seeds 1-10, dropout 0/0.1/0.3 nested masks seed 5000+rep).
# Part B: real UCSF CD8T pseudobulk (99 donors, top-2000 variable genes,
#         log1p-CP10K from .raw), 3 seeds.
#
# Arms: S=NULL (package default) vs S=vignette floor (n = nrow of the fit
# matrix). Metrics: K retained, score recovery (Hungarian |cor|), top-30
# Jaccard, lfsr AUC (planted vs noise genes), false lfsr<0.05 rate among
# noise genes, mean estimated residual SD vs S, ELBO.
#
# Audit pilot — NOT the registered MR-2 benchmark. No R/ package change here.

set.seed(20261018)
suppressPackageStartupMessages({
  library(flashier); library(ebnm)
})
options(cellprograms.use_flashier_utils = FALSE)  # direct flashier: S + lfsr access

SEED <- 20261018
outdir <- "/workspace/cand13"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- generator (E-workstream recovery_k3, pinned) --------------------------
gen_k3 <- function(seed, n = 80, p = 600, block = 30, effect = 1.0) {
  set.seed(seed)
  K <- 3
  z <- matrix(rnorm(n * K), n, K)
  z <- z + matrix(rnorm(n * K, sd = 0.3), n, K)   # donor noise sd 0.3
  W <- matrix(0, p, K)
  for (k in seq_len(K)) W[((k - 1) * block + 1):(k * block), k] <- effect
  Y <- z %*% t(W) + matrix(rnorm(n * p), n, p)
  list(Y = Y, z = z, W = W)
}
gen_dropout_mask <- function(seed, Y, rate) {
  if (rate <= 0) return(matrix(TRUE, nrow(Y), ncol(Y)))
  set.seed(seed)
  u <- matrix(runif(length(Y)), nrow(Y))
  mask <- u > rate
  # zeros replace the masked entries on the COUNT scale; the latent Y is
  # log-scale, so the mask is applied to counts = round(exp(Y)) downstream.
  mask
}

shifted_log <- function(counts) {
  sf <- rowSums(counts); sf <- sf / mean(sf)
  log1p(sweep(counts, 1, sf, "/"))
}

hungarian_match <- function(est, truth) {
  # max |cor| bipartite match via greedy on the cor matrix (K<=10, exact
  # assignment not needed for calibration metrics; E used the same greedy)
  C <- abs(cor(est, truth))
  pairs <- list(); used_e <- used_t <- c()
  while (length(used_e) < min(ncol(est), ncol(truth))) {
    mx <- which(C == max(C), arr.ind = TRUE)[1, ]
    if (C[mx[1], mx[2]] <= 0) break
    pairs[[length(pairs) + 1]] <- c(est_k = mx[1], true_k = mx[2], cor = C[mx[1], mx[2]])
    used_e <- c(used_e, mx[1]); used_t <- c(used_t, mx[2])
    C[mx[1], ] <- -1; C[, mx[2]] <- -1
  }
  do.call(rbind, pairs)
}

vignette_S <- function(n) {
  x <- rpois(1e7, 1 / n)
  sd(log(x + 1))
}

fit_arm <- function(Ylog, S, max_factors = 10, seed = 1) {
  set.seed(seed)
  args <- list(data = Ylog, ebnm_fn = list(ebnm::ebnm_normal, ebnm::ebnm_point_laplace),
               var_type = 1, greedy_Kmax = max_factors, backfit = TRUE,
               nullcheck = TRUE, verbose = 0)
  if (!is.null(S)) args$S <- S
  fl <- tryCatch(do.call(flashier::flash, args), error = function(e) NULL)
  if (is.null(fl)) return(list(K = NA, fit = NULL))
  list(K = fl$n_factors, fit = fl)
}

lfsr_stats <- function(fl, W_true, block = 30, K = 3) {
  # AUC of 1 - lfsr separating planted from noise genes (pooled over matched
  # factors); false rate = fraction of noise genes with lfsr < 0.05.
  Fl <- fl$F_lfsr
  if (is.null(Fl) || fl$n_factors == 0) return(c(auc = NA, false_rate = NA))
  planted <- rep(FALSE, nrow(Fl))
  for (k in seq_len(K)) planted[((k - 1) * block + 1):(k * block)] <- TRUE
  # match factors to truth blocks by |cor| of loadings
  C <- abs(cor(fl$F_pm, W_true))
  matched <- if (is.null(dim(C))) c() else apply(C, 1, which.max)
  matched <- unique(matched[matched <= K])
  if (!length(matched)) return(c(auc = NA, false_rate = NA))
  lf <- as.vector(Fl[, matched, drop = FALSE])
  lab <- rep(planted, length(matched))
  lf1 <- 1 - lf
  # AUC via rank statistic
  r <- rank(lf1)
  n1 <- sum(lab); n0 <- sum(!lab)
  auc <- (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  c(auc = auc, false_rate = mean(lf[!lab] < 0.05, na.rm = TRUE))
}

# ---- Part A: dropout scenarios ---------------------------------------------
rows <- list()
for (rate in c(0, 0.1, 0.3)) {
  S_v <- vignette_S(80)
  for (rep in 1:10) {
    g <- gen_k3(seed = rep)
    counts <- round(exp(g$Y))
    mask <- gen_dropout_mask(seed = 5000 + rep, g$Y, rate)
    counts_d <- counts * mask
    Ylog <- shifted_log(counts_d)
    for (arm in c("default", "vignetteS")) {
      S <- if (arm == "vignetteS") S_v else NULL
      fa <- fit_arm(Ylog, S, seed = rep)
      if (is.null(fa$fit)) {
        rows[[length(rows) + 1]] <- data.frame(part = "A", rate = rate, rep = rep,
          arm = arm, K = NA, score_cor = NA, jaccard = NA, auc = NA,
          false_rate = NA, mean_sd = NA, S = S_v, elbo = NA)
        next
      }
      fl <- fa$fit
      sc <- NA
      if (fl$n_factors > 0) {
        m <- hungarian_match(fl$L_pm, g$z)
        sc <- if (is.null(m)) NA else max(m[, "cor"])
      }
      # top-30 Jaccard pooled over matched factors
      jac <- NA
      if (fl$n_factors > 0) {
        C <- abs(cor(fl$F_pm, g$W))
        matched <- unique(apply(C, 1, which.max))
        matched <- matched[matched <= 3]
        if (length(matched)) {
          js <- sapply(matched, function(k) {
            top <- order(-abs(fl$F_pm[, k]))[1:30]
            truth <- which(g$W[, k] != 0)
            length(intersect(top, truth)) / length(union(top, truth))
          })
          jac <- mean(js)
        }
      }
      ls <- if (fl$n_factors > 0) lfsr_stats(fl, g$W) else c(auc = NA, false_rate = NA)
      rows[[length(rows) + 1]] <- data.frame(part = "A", rate = rate, rep = rep,
        arm = arm, K = fa$K, score_cor = sc, jaccard = jac,
        auc = ls[["auc"]], false_rate = ls[["false_rate"]],
        mean_sd = mean(fl$residuals_sd), S = S_v,
        elbo = fl$elbo)
    }
  }
  cat("Part A rate", rate, "done\n"); flush.console()
  write.csv(do.call(rbind, rows), file.path(outdir, "cand13_partA_partial.csv"), row.names = FALSE)
}

# ---- Part B: real UCSF CD8T pseudobulk -------------------------------------
# counts from .raw (raw counts layer), CD8 T cells, 99 donors.
suppressPackageStartupMessages({ library(Matrix) })
h5 <- "/workspace/data/UCSF_reannotated.h5ad"
if (!file.exists(h5)) h5 <- "/mnt/shared-workspace/shared/audit/data/UCSF_reannotated.h5ad"
has_h5 <- requireNamespace("rhdf5", quietly = TRUE)
if (has_h5) {
  cat("Part B: reading", h5, "\n"); flush.console()
  # anndata raw/X is CSR (encoding-type csr_matrix): rows = cells; pieces at
  # raw/X/{data,indices,indptr}; shape is an ATTRIBUTE of the raw/X group.
  xd <- rhdf5::h5read(h5, "raw/X/data")
  xi <- rhdf5::h5read(h5, "raw/X/indices")
  xp <- rhdf5::h5read(h5, "raw/X/indptr")
  xsh <- as.integer(rhdf5::h5readAttributes(h5, "raw/X")$shape)  # c(n_cells, n_genes)
  X <- new("dgRMatrix", j = as.integer(xi), p = as.integer(xp),
           x = as.numeric(xd), Dim = xsh)
  rm(xd, xi, xp); gc()
  # categoricals: {categories, codes} pairs (newer anndata) or raw datasets
  read_cat <- function(path) {
    g <- rhdf5::h5read(h5, path)
    if (is.list(g) && !is.null(g$categories)) {
      cats <- as.character(g$categories)
      codes <- as.vector(g$codes)
      out <- cats[codes + 1L]
      out[codes < 0] <- NA
      out
    } else as.character(g)
  }
  cd8 <- read_cat("obs/author_celltype")
  donors <- read_cat("obs/donor_id")
  genes_all <- as.character(rhdf5::h5read(h5, "raw/var/_index"))
  keep <- which(!is.na(cd8) & cd8 == "CD8 T cell" & !is.na(donors))
  d_kept <- donors[keep]
  ud <- sort(unique(d_kept))
  cat("CD8 T cells:", length(keep), "donors:", length(ud), "\n"); flush.console()
  Xc <- X[keep, , drop = FALSE]
  colnames(Xc) <- genes_all
  rm(X); gc()
  d_f <- factor(d_kept, levels = ud)
  G <- Matrix::Matrix(0, nrow = length(d_f), ncol = length(ud), sparse = TRUE)
  G[cbind(seq_along(d_f), as.integer(d_f))] <- 1
  pb <- as.matrix(Matrix::t(G) %*% Matrix::t(Xc))   # donors x genes (raw counts)
  colnames(pb) <- colnames(Xc); rownames(pb) <- ud
  # top-2000 variable genes on log1p-CP10K
  lib <- rowSums(pb)
  cp10k <- log1p(sweep(pb, 1, lib, function(x, s) x / s * 1e4))
  vv <- apply(cp10k, 2, var)
  genes <- names(sort(vv, decreasing = TRUE))[1:2000]
  Yb <- cp10k[, genes]
  rm(Xc, cp10k); gc()
  S_v <- vignette_S(nrow(Yb))
  for (seed in c(20261016, 20261017, 20261018)) {
    for (arm in c("default", "vignetteS")) {
      S <- if (arm == "vignetteS") S_v else NULL
      t0 <- Sys.time()
      fa <- fit_arm(Yb, S, max_factors = 10, seed = seed)
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(fa$fit)) {
        rows[[length(rows) + 1]] <- data.frame(part = "B", rate = NA, rep = seed,
          arm = arm, K = NA, score_cor = NA, jaccard = NA, auc = NA,
          false_rate = NA, mean_sd = NA, S = S_v, elbo = NA)
        next
      }
      fl <- fa$fit
      Fl <- fl$F_lfsr
      fr <- if (is.null(Fl) || fl$n_factors == 0) NA else mean(as.vector(Fl) < 0.05, na.rm = TRUE)
      rows[[length(rows) + 1]] <- data.frame(part = "B", rate = NA, rep = seed,
        arm = arm, K = fa$K, score_cor = NA, jaccard = NA, auc = NA,
        false_rate = fr, mean_sd = mean(fl$residuals_sd), S = S_v,
        elbo = fl$elbo)
      # save loadings for seed-stability across seeds
      if (fl$n_factors > 0) write.csv(fl$F_pm, file.path(outdir,
        sprintf("B_%s_%d_loadings.csv", arm, seed)), row.names = genes)
      cat("Part B", arm, seed, "K =", fa$K, "in", round(el, 1), "s\n"); flush.console()
    }
  }
} else {
  cat("Part B SKIPPED: rhdf5 not available\n")
}

res <- do.call(rbind, rows)
write.csv(res, file.path(outdir, "cand13_demonstrator.csv"), row.names = FALSE)
cat("DONE. rows:", nrow(res), "\n")
cat("S vignette (n=80):", vignette_S(80), " S vignette (n=99):", vignette_S(99), "\n")
