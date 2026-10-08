## 04_corrected_eval.R <cohort> <fit_rds|label=rds[,..]|NONE> <out_csv>
## Corrected metrics, independent of Biomni's:
##  (1) dimension-matched random floor: Gaussian rep with the SAME number of columns, mean over 20 draws.
##      `floor`        = fully observed Gaussian.
##      `floor_masked` = Gaussian with the representation's own missing-donor pattern (zeros / blocks masked),
##                       i.e. the floor that also absorbs missingness-pattern leakage.
##  (2) available-case kNN: squared distance summed only over blocks observed in BOTH donors, rescaled by
##      B / (#shared blocks); pairs with < 2 shared blocks are dropped (single-block reps: min = 1).
##  (3) 10-fold CV ridge R^2 (lambda by exact LOO inside each training fold, via SVD) replaces max|Spearman|.
## Output: representation,covariate,type,metric,value,panel,floor,floor_masked,metric_variant,n_dim

source(file.path(local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
                         if (is.na(f)) getwd() else dirname(normalizePath(f)) }), "02_metrics.R"))

## ---- (2) available-case distance (squared, B/shared rescale); Inf = dropped pair
ac_dist <- function(X, blocks, observed, min_shared = 2L) {
  B <- ncol(observed); ms <- min(min_shared, B)
  X <- unclass_matrix(X)
  Dsum <- 0; cnt <- 0
  for (b in colnames(observed)) {
    Xb <- X[, blocks == b, drop = FALSE]
    r2 <- rowSums(Xb^2)
    d2 <- pmax(outer(r2, r2, "+") - 2 * tcrossprod(Xb), 0)
    M <- tcrossprod(as.numeric(observed[, b]))
    Dsum <- Dsum + d2 * M; cnt <- cnt + M     # unobserved rows are masked by M, whatever they hold
  }
  D <- Dsum / cnt * B
  D[cnt < ms] <- Inf; diag(D) <- Inf
  D
}

knn_ac_macro_f1 <- function(D, y, k = 5L) {
  ok <- which(!is.na(y)); if (length(ok) < k + 2 || length(unique(y[ok])) < 2) return(NA_real_)
  D <- D[ok, ok, drop = FALSE]; yy <- as.character(y[ok])
  pred <- vapply(seq_along(ok), function(i) {
    o <- order(D[i, ]); o <- o[is.finite(D[i, o])][seq_len(k)]; o <- o[!is.na(o)]
    if (!length(o)) return(names(which.max(table(yy[-i]))))
    v <- table(yy[o]); w <- names(v)[v == max(v)]
    if (length(w) == 1L) w else yy[o][yy[o] %in% w][1]   # tie -> class of the nearest tied neighbour
  }, character(1))
  macro_f1(yy, pred)
}

## ---- (3) CV ridge R^2.  Intercept via train-fold centring; the LOO shortcut is exact for that model only if the
## leverage includes the intercept's 1/n_tr term (omitting it made h->1 never bind, so lambda->min won when
## p ~ n_tr, the fit interpolated and out-of-fold R^2 was hugely negative).
cv_ridge_r2 <- function(X, y, nfold = 10L, seed = 1L) {
  ok <- !is.na(y); X <- unclass_matrix(X)[ok, , drop = FALSE]; y <- y[ok]; n <- length(y)
  if (n < nfold * 2 || stats::sd(y) == 0) return(NA_real_)
  set.seed(seed); fold <- sample(rep_len(seq_len(nfold), n)); pred <- numeric(n)
  for (f in seq_len(nfold)) {
    tr <- fold != f; mu <- colMeans(X[tr, , drop = FALSE]); yb <- mean(y[tr])
    Xc <- sweep(X[tr, , drop = FALSE], 2, mu); yc <- y[tr] - yb
    s <- svd(Xc); d2 <- s$d^2
    if (max(d2) < 1e-10) { pred[!tr] <- yb; next }
    uy <- crossprod(s$u, yc); u2 <- s$u^2
    ## penalty floor 1e-4*max(d2): when p >= n_tr the tiny-lambda end interpolates the training fold (LOO is a poor
    ## selector there: n small, noisy, near-singular); a floor keeps a minimum of shrinkage as a safeguard.
    lam <- max(d2) * 10^seq(-4, 1, length.out = 30)
    loo <- vapply(lam, function(l) {
      w <- d2 / (d2 + l); h <- 1 / sum(tr) + as.numeric(u2 %*% w)   # + intercept leverage 1/n_tr
      mean(((yc - s$u %*% (w * uy)) / (1 - h))^2)
    }, numeric(1))
    l <- lam[which.min(loo)]
    beta <- s$v %*% ((s$d / (d2 + l)) * uy)
    pred[!tr] <- yb + sweep(X[!tr, , drop = FALSE], 2, mu) %*% beta
  }
  1 - sum((y - pred)^2) / sum((y - mean(y))^2)
}

## ---- all corrected + reference metrics for ONE (blocked) representation; returns long data.frame
eval_rep_corrected <- function(R, meta, spec, k = 5L, seed = 1L) {
  set.seed(seed)
  blocks <- attr(R, "blocks"); obs <- attr(R, "observed")
  D <- ac_dist(R, blocks, obs)
  cat_cols <- intersect(c(spec$cat_bio, spec$cat_tech), names(meta))
  cont_cols <- intersect(spec$cont_bio, names(meta))
  rows <- list()
  for (cv in cat_cols) {
    rows[[length(rows) + 1L]] <- data.frame(covariate = cv, type = "categorical", metric = "knn_macro_f1",
      metric_variant = c("biomni_zerofill", "availcase"),
      value = c(knn_macro_f1(R, meta[[cv]], k), knn_ac_macro_f1(D, meta[[cv]], k)))
  }
  for (cv in cont_cols) {
    y <- as.numeric(meta[[cv]])
    rows[[length(rows) + 1L]] <- data.frame(covariate = cv, type = "continuous",
      metric = c("max_abs_spearman", "cv_ridge_r2"), metric_variant = c("biomni_spearman", "cv_ridge"),
      value = c(max_abs_spearman(R, y), cv_ridge_r2(R, y)))
  }
  do.call(rbind, rows)
}

random_like <- function(R, masked) {
  X <- matrix(stats::rnorm(nrow(R) * ncol(R)), nrow(R), ncol(R), dimnames = dimnames(R))
  obs <- attr(R, "observed"); blocks <- attr(R, "blocks")
  if (masked) for (b in colnames(obs)) X[!obs[, b], blocks == b] <- 0
  else obs[] <- TRUE
  attr(X, "blocks") <- blocks; attr(X, "observed") <- obs
  X
}

floor_values <- function(R, meta, spec, ndraw = 20L, masked = FALSE, seed = 1000L) {
  draws <- lapply(seq_len(ndraw), function(i) {
    set.seed(seed + i)
    eval_rep_corrected(random_like(R, masked), meta, spec, seed = seed + i)$value
  })
  rowMeans(do.call(cbind, draws), na.rm = TRUE)
}

eval_corrected <- function(reps, spec, donors, ndraw = as.integer(Sys.getenv("BIOMNI_NDRAW", "20"))) {
  meta <- spec$meta[donors, , drop = FALSE]; rownames(meta) <- donors
  do.call(rbind, lapply(names(reps), function(rn) {
    R <- reps[[rn]]
    ev <- eval_rep_corrected(R, meta, spec)
    ev$floor <- floor_values(R, meta, spec, ndraw, masked = FALSE)
    ev$floor_masked <- floor_values(R, meta, spec, ndraw, masked = TRUE)
    data.frame(representation = rn, ev, n_dim = ncol(R))
  }))
}

finish_eval <- function(res, spec) {
  res$panel <- ifelse(res$covariate %in% spec$cat_tech, "B_technical", "A_biology")
  res <- res[, c("representation", "covariate", "type", "metric", "value", "panel", "floor", "floor_masked",
                 "metric_variant", "n_dim")]
  rownames(res) <- NULL; res
}

if (!interactive() && identical(basename(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                                "04_corrected_eval.R")) {
  a <- commandArgs(TRUE)
  if (length(a) != 3) stop("usage: Rscript 04_corrected_eval.R <cohort> <fit_rds|NONE> <out_csv>")
  fits <- parse_fits(a[2])
  if (length(fits)) pkgload::load_all(file.path(.script_dir(), "fork"), quiet = TRUE)
  L <- load_cohort(a[1]); spec <- cohort_spec(L$meta, a[1])
  res <- finish_eval(eval_corrected(all_representations(L, fits), spec, L$donors), spec)
  utils::write.csv(res, a[3], row.names = FALSE)
  cat("wrote", nrow(res), "rows to", a[3], "\n")
}
