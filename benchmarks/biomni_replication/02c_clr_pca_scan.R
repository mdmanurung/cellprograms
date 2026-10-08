## 02c_clr_pca_scan.R -- stage 2 of the variant search: clr pseudocount and PCA preprocessing, scored on the cells
## that are (near-)deterministic (binary covariates: no kNN vote ties) plus the continuous Spearman rows.
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
source(file.path(here, "02b_metrics_variants.R"))
donors <- L$donors; meta <- sp$meta[donors, ]
refv <- setNames(ref$value, paste(ref$representation, ref$covariate))
det <- c("Death28", "sex", "Institute"); cont <- sp$cont_bio
V0 <- list(dist = "euclid", tie = "lowclass", use_all = FALSE, l = 0L, na = "drop")
score <- function(X, rn, covs = c(det, "Source", "Outcome", "GEX_region")) {
  k <- sapply(covs, function(cv) cell_value(X, meta[[cv]], V0)["truth"]); names(k) <- covs
  cs <- sapply(cont, function(cv) max_abs_spearman(X, as.numeric(meta[[cv]]))); names(cs) <- cont
  r <- refv[paste(rn, c(covs, cont))]
  d <- c(k, cs) - r
  c(det_mad = mean(abs(d[det])), all_cat_mad = mean(abs(d[covs])), cont_mad = mean(abs(d[cont])), round(d, 3))
}
options(width = 220)
p <- L$comp / rowSums(L$comp); p <- p[donors, ]
cat("== clr (unscaled), pseudocount scan; ref = clr_composition\n")
res <- NULL
for (pc in c(1e-8, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 0.1, 0.5, 1))
  res <- rbind(res, data.frame(variant = paste0("p+", pc), t(score({l <- log(p + pc); l - rowMeans(l)}, "clr_composition"))))
for (f in c(0.5, 1)) { m <- min(p[p > 0]); q <- p; q[q == 0] <- f * m; l <- log(q / rowSums(q)); res <- rbind(res, data.frame(variant = paste0("zero->", f, "*min"), t(score(l - rowMeans(l), "clr_composition")))) }
for (n in c(100, 1000, 5000, 10000)) { q <- p * n + 0.5; q <- q / rowSums(q); l <- log(q); res <- rbind(res, data.frame(variant = paste0("counts~", n, "+0.5"), t(score(l - rowMeans(l), "clr_composition")))) }
l <- log(p + 1e-4); l <- l - rowMeans(l); res <- rbind(res, data.frame(variant = "p+1e-4, z-scored", t(score(zs(l), "clr_composition"))))
print(res[order(res$det_mad), ], digits = 3, row.names = FALSE)

cat("\n== PCA preprocessing, global_pca (ref) -- fully observed\n")
pca2 <- function(X, k, scale = FALSE, log = FALSE, center = TRUE) { X <- as.matrix(X); if (log) X <- log1p(X)
  X <- X[, apply(X, 2, sd) > 1e-12, drop = FALSE]; prcomp(X, center = center, scale. = scale)$x[, seq_len(k), drop = FALSE] }
res <- NULL
for (sc in c(FALSE, TRUE)) for (k in c(5, 10, 20, 30)) for (zz in c("none", "z")) {
  P <- pca2(L$global, k, sc); if (zz == "z") P <- zs(P)
  res <- rbind(res, data.frame(variant = sprintf("prcomp scale=%s k=%d pcs=%s", sc, k, zz), t(score(P[donors, ], "global_pca"))))
}
print(res[order(res$det_mad), ][1:10, ], digits = 3, row.names = FALSE)
