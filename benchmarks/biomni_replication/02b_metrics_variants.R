## 02b_metrics_variants.R -- diagnostic grid: which defensible metrics.R variant reproduces Biomni's COMBAT
## categorical 5-NN macro-F1 for the baselines?  Usage: BIOMNI_DATA=<data> Rscript 02b_metrics_variants.R
## Writes results/diag/knn_variants_all.csv (every variant x rep x covariate) and knn_variants.csv (top 5).
## Does not modify 02_metrics.R (sourced for load_cohort / cohort_spec / .safe_pca / rep_random).

here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
source(file.path(here, "02_metrics.R"))
dir.create(file.path(here, "results/diag"), showWarnings = FALSE, recursive = TRUE)
NSIM <- 20L  # replicates for "random" tie-breaking (expectation of macro-F1)

L <- load_cohort("combat"); sp <- cohort_spec(L$meta)
ref <- read.csv(file.path(here, "ref/cellprograms_benchmark/combat/tables/eval_all.csv"))
cats <- sp$cat_bio; cats <- c(cats, sp$cat_tech)
reps_scored <- c("composition", "clr_composition", "global_pca", "perct_pca", "grouped_pca")

## ----------------------------------------------------------------- representation builders (raw, observed-only)
zs <- function(M) { M <- as.matrix(M); s <- apply(M, 2, sd); s[!is.finite(s) | s < 1e-12] <- 1
  sweep(sweep(M, 2, colMeans(M), "-"), 2, s, "/") }
fill0 <- function(sc, donors) { o <- matrix(0, length(donors), ncol(sc), dimnames = list(donors, colnames(sc)))
  cm <- intersect(rownames(sc), donors); o[cm, ] <- sc[cm, , drop = FALSE]; o }
## scale = none | zbefore (z on observed donors, then zero-fill) | zafter (zero-fill, then z over all donors)
assemble <- function(blocks, donors, scale) {
  do.call(cbind, lapply(blocks, function(m) switch(scale,
    none = fill0(m, donors), zbefore = fill0(zs(m[rownames(m) %in% donors, , drop = FALSE]), donors),
    zafter = zs(fill0(m, donors)))))
}
make_reps <- function(donors, comp_scale, pca_scale) {
  p <- L$comp / rowSums(L$comp)
  cl <- { l <- log(p + 1e-4); l - rowMeans(l) }
  cs <- function(M) switch(comp_scale, none = M[donors, , drop = FALSE], z = zs(M)[donors, , drop = FALSE],
                           z_obs = zs(M[donors, , drop = FALSE]))
  list(composition = cs(p), clr_composition = cs(cl),
       global_pca = assemble(list(.safe_pca(L$global, 10L)), donors, pca_scale),
       perct_pca = assemble(lapply(L$mats, .safe_pca, k = 3L), donors, pca_scale),
       grouped_pca = assemble(lapply(L$groups, .safe_pca, k = 5L), donors, pca_scale))
}

## ----------------------------------------------------------------- kNN (manual, so tie rules are explicit)
dmat <- function(X, dist) {
  if (dist == "cosine") { n <- sqrt(rowSums(X^2)); n[n < 1e-12] <- 1; Xn <- X / n; D <- 1 - Xn %*% t(Xn); D[D < 0] <- 0 }
  else { s <- rowSums(X^2); D <- outer(s, s, "+") - 2 * X %*% t(X); D[D < 0] <- 0; D <- sqrt(D) }
  diag(D) <- Inf; D
}
## tie: random | nearest (tied class of the closest neighbour) | lowclass (first level, which.max) ; use_all: keep
## neighbours tied with the k-th distance; l: minimum winning votes else NA ("doubt")
knn_pred <- function(D, O, y, k = 5L, tie = "random", use_all = FALSE, l = 0L) {
  n <- nrow(D); lev <- levels(y); yi <- as.integer(y); pred <- integer(n)
  for (i in seq_len(n)) {
    o <- O[i, ]; d <- D[i, o]
    nb <- if (use_all) o[is.finite(d) & d <= d[k] * (1 + 1e-9) + 1e-12] else o[seq_len(k)]
    v <- tabulate(yi[nb], length(lev)); w <- which(v == max(v))
    pred[i] <- if (length(w) == 1L) w else switch(tie, lowclass = w[1],
      nearest = { cand <- yi[nb][yi[nb] %in% w]; cand[1] },
      random = w[sample.int(length(w), 1L)])
    if (v[pred[i]] < l) pred[i] <- NA_integer_
  }
  factor(lev[pred], levels = lev)
}

## macro-F1 variants
f1_one <- function(t, p, cl) { d <- sum(t == cl) + sum(p == cl, na.rm = TRUE); if (d == 0) 0 else 2 * sum(t == cl & p == cl, na.rm = TRUE) / d }
macro_f1v <- function(t, p, def) {
  t <- as.character(t); p <- as.character(p)
  switch(def,
    truth = mean(vapply(unique(t), function(c) f1_one(t, p, c), 0)),                          # current (fixed Biomni)
    union = mean(vapply(union(t, na.omit(p)), function(c) f1_one(t, p, c), 0)),
    pred_only = mean(vapply(unique(na.omit(p)), function(c) f1_one(t, p, c), 0)),             # pre-fix style: classes predicted
    avgPR = { cl <- unique(t); pr <- vapply(cl, function(c) { n <- sum(p == c, na.rm = TRUE); if (n == 0) 0 else sum(t == c & p == c, na.rm = TRUE) / n }, 0)
      rc <- vapply(cl, function(c) sum(t == c & p == c, na.rm = TRUE) / sum(t == c), 0)
      P <- mean(pr); R <- mean(rc); if (P + R == 0) 0 else 2 * P * R / (P + R) })
}

## one cell: expected macro-F1 (all F1 definitions at once) for (rep matrix, labels) under a kNN variant
F1DEFS <- c("truth", "avgPR", "union", "pred_only")
cell_value <- function(X, y, v) {
  y <- as.character(y)
  if (v$na == "keep") y[is.na(y)] <- "NA_label" else { ok <- !is.na(y); X <- X[ok, , drop = FALSE]; y <- y[ok] }
  if (length(unique(y)) < 2) return(setNames(rep(NA_real_, length(F1DEFS)), F1DEFS))
  yf <- factor(y); D <- dmat(X, v$dist); O <- t(apply(D, 1, order))
  nsim <- if (v$tie == "random") NSIM else 1L
  rowMeans(vapply(seq_len(nsim), function(s) {
    p <- knn_pred(D, O, yf, 5L, v$tie, v$use_all, v$l)
    vapply(F1DEFS, function(d) macro_f1v(y, p, d), 0)
  }, numeric(length(F1DEFS))))
}

if (!nzchar(Sys.getenv("VARIANTS_LIB"))) {  # set VARIANTS_LIB=1 to source only the functions
## ----------------------------------------------------------------- variant grid
grid <- expand.grid(comp_scale = c("none", "z", "z_obs"), pca_scale = c("none", "zbefore", "zafter"),
                    dist = c("euclid", "cosine"),
                    knn = c("random", "nearest", "lowclass", "random_useall", "random_l3"),
                    na = c("drop", "keep"),
                    donors = c("122", "124"), stringsAsFactors = FALSE)
## z_obs == z when donors == 124 -> drop duplicates; NA handling only matters for Outcome (kept cheap anyway)
grid <- grid[!(grid$comp_scale == "z_obs" & grid$donors == "124"), ]
grid$knn_tie <- ifelse(grid$knn %in% c("random", "random_useall", "random_l3"), "random", grid$knn)
grid$use_all <- grid$knn == "random_useall"; grid$l <- ifelse(grid$knn == "random_l3", 3L, 0L)
cat("variants:", nrow(grid), "\n")

meta124 <- sp$meta[sort(rownames(L$comp)), , drop = FALSE]
refv <- setNames(ref$value, paste(ref$representation, ref$covariate))
set.seed(1)
cache <- list()
rows <- vector("list", nrow(grid))
for (g in seq_len(nrow(grid))) {
  v <- grid[g, ]
  donors <- if (v$donors == "122") L$donors else sort(union(L$donors, rownames(L$comp)))
  key <- paste(v$comp_scale, v$pca_scale, v$donors)
  if (is.null(cache[[key]])) cache[[key]] <- make_reps(donors, v$comp_scale, v$pca_scale)
  reps <- cache[[key]]
  meta <- sp$meta[donors, , drop = FALSE]
  vv <- list(dist = v$dist, tie = v$knn_tie, use_all = v$use_all, l = v$l, na = v$na)
  out <- list()
  for (rn in reps_scored) for (cv in cats)
    out[[length(out) + 1L]] <- data.frame(vid = g, representation = rn, covariate = cv, f1 = F1DEFS,
      value = cell_value(reps[[rn]], meta[[cv]], vv), ref = refv[[paste(rn, cv)]])
  rows[[g]] <- do.call(rbind, out)
  if (g %% 20 == 0) cat(g, "")
}
cat("\n")
all <- do.call(rbind, rows)
all$absdiff <- abs(all$value - all$ref)
grid$vid <- seq_len(nrow(grid))
sc <- aggregate(absdiff ~ vid + f1, all, mean); names(sc)[3] <- "mean_abs_diff"
n_ok <- aggregate(absdiff ~ vid + f1, all, function(x) sum(x < 0.03)); names(n_ok)[3] <- "n_cells_lt_0.03"
summ <- merge(merge(grid, sc), n_ok)
summ$n_departures <- with(summ, (comp_scale != "none") + (pca_scale != "none") + (dist != "euclid") + (knn != "random") +
                           (f1 != "truth") + (na != "drop") + (donors != "122"))
summ <- summ[order(summ$mean_abs_diff, summ$n_departures), ]
write.csv(all, file.path(here, "results/diag/knn_variants_all.csv"), row.names = FALSE)
write.csv(summ, file.path(here, "results/diag/knn_variants_summary.csv"), row.names = FALSE)
cur <- subset(summ, comp_scale == "z" & pca_scale == "zbefore" & dist == "euclid" & knn == "random" & f1 == "truth" & na == "drop" & donors == "122")
cat("current implementation:\n"); print(cur[, c("mean_abs_diff", "n_cells_lt_0.03")])
cat("top 15:\n"); print(head(summ[, c("vid", "comp_scale", "pca_scale", "dist", "knn", "f1", "na", "donors", "mean_abs_diff", "n_cells_lt_0.03", "n_departures")], 15), row.names = FALSE)
saveRDS(list(grid = grid, summ = summ), file.path(here, "results/diag/knn_variants.rds"))
}
