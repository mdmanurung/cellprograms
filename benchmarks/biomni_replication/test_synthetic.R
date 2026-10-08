## test_synthetic.R -- unit checks of 02/03/04 on a synthetic COMBAT-shaped cohort (needs the cellprograms-r env).
## Usage: Rscript test_synthetic.R   (writes results/synth_data/combat_synth, results/synth_out/*)
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
Sys.setenv(BIOMNI_DATA = file.path(here, "results/synth_data"), BIOMNI_NDRAW = "3")
RS <- file.path(R.home("bin"), "Rscript"); out <- file.path(here, "results/synth_out"); dir.create(out, TRUE, TRUE)
run <- function(...) { r <- system2(RS, c(...), stdout = TRUE, stderr = TRUE); st <- attr(r, "status")
  if (!is.null(st) && st != 0) { cat(r, sep = "\n"); stop("command failed: ", paste(...)) }; invisible(r) }

## ---- synthetic cohort ----
set.seed(7); N <- 70; donors <- sprintf("D%03d", 1:N); G <- 60
d <- file.path(here, "results/synth_data/combat_synth"); dir.create(file.path(d, "mat"), TRUE, TRUE)
sex <- sample(c("male", "female"), N, TRUE); age <- runif(N, 20, 90)
band <- cut(age, c(0, 30, 40, 50, 60, 70, 80, 90, 200), labels = c("19-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", ">=91"))
inst <- sample(c("Oxford", "StG"), N, TRUE, prob = c(.85, .15))
meta <- data.frame(Source = sample(c("HV", "SEV", "MILD"), N, TRUE), Outcome = sample(c("a", "b", "c", "d"), N, TRUE),
                   Death28 = sample(c("True", "False"), N, TRUE, prob = c(.2, .8)), sex = sex, Institute = inst,
                   GEX_region = sample(c("r1", "r2", "r3"), N, TRUE), Age = as.character(band),
                   Hospitalstay = round(rnorm(N, 8, 4)), TimeSinceOnset = round(rnorm(N, 10, 3)), row.names = donors)
meta$Hospitalstay[sample(N, 15)] <- NA; meta$Age[sample(N, 3)] <- NA
cts <- c("CD4", "B", "NK", "Mono"); mats <- list()
for (ct in cts) {
  keep <- sort(sample(donors, sample(45:62, 1)))
  Z <- cbind(ifelse(sex[match(keep, donors)] == "male", 1, -1), scale(age[match(keep, donors)])[, 1])
  W <- matrix(rnorm(2 * G), 2) * (runif(2 * G) < .3) * 1.5
  M <- Z %*% W + matrix(rnorm(length(keep) * G, sd = .8), length(keep)); dimnames(M) <- list(keep, paste0("g", 1:G))
  mats[[ct]] <- M; write.csv(M, file.path(d, "mat", paste0("ct_", ct, ".csv")))
}
allM <- do.call(rbind, lapply(mats, function(m) m[, 1:20])); glob <- rowsum(allM, rownames(allM)) / as.vector(table(rownames(allM)))
write.csv(glob, file.path(d, "mat", "global.csv"))
ci <- intersect(rownames(mats$CD4), rownames(mats$NK)); write.csv(cbind(mats$CD4[ci, ], mats$NK[ci, ]), file.path(d, "mat", "group_T.csv"))
write.csv(mats$Mono, file.path(d, "mat", "group_Mono.csv"))
comp <- matrix(rpois(N * 6, 30), N, 6, dimnames = list(donors, paste0("ct", 1:6))); write.csv(comp, file.path(d, "composition.csv"))
write.csv(meta, file.path(d, "donor_meta.csv"))

## ---- primitive unit checks (sourced libs) ----
source(file.path(here, "04_corrected_eval.R"))
stopifnot(abs(macro_f1(c(rep("O", 108), rep("S", 14)), rep("O", 122)) - 0.4698) < 1e-3)  # constant predictor -> ~0.47 floor
X <- matrix(rnorm(80 * 12), 80, dimnames = list(sprintf("d%d", 1:80), NULL)); y <- ifelse(X[, 1] > 0, "u", "v")
R <- blocked_rep(list(a = X[, 1:4], b = X[, 5:8], c = X[, 9:12]), rownames(X))
D <- ac_dist(R, attr(R, "blocks"), attr(R, "observed"))
stopifnot(max(abs(D[upper.tri(D)] - as.matrix(dist(unclass_matrix(R)))[upper.tri(D)]^2)) < 1e-8)  # fully observed == Euclid^2
set.seed(1); stopifnot(abs(knn_ac_macro_f1(D, y) - knn_macro_f1(R, y)) < 0.1)
Xm <- R; ob <- attr(R, "observed"); ob[1:10, "a"] <- FALSE; ob[1:10, "b"] <- FALSE   # donors 1-10 have 1 block -> dropped pairs
Xm[1:10, 1:8] <- 0; attr(Xm, "observed") <- ob
Dm <- ac_dist(Xm, attr(R, "blocks"), ob); stopifnot(all(is.infinite(Dm[1:10, 1:10])), all(is.finite(Dm[11:20, 11:20][upper.tri(diag(10))])))
sig <- X[, 1:3] %*% c(1, 1, 1) + rnorm(80, sd = .5); noise <- rnorm(80)
r1 <- cv_ridge_r2(X, sig); r0 <- cv_ridge_r2(X, noise); cat(sprintf("ridge R2 signal %.2f, noise %.2f\n", r1, r0))
stopifnot(r1 > 0.6, r0 < 0.1)

## ---- regression: ridge LOO leverage must include the intercept (p ~ n_train used to interpolate, R2 ~ -20) ----
## OLD (buggy) version kept ONLY to demonstrate the failure: leverage without 1/n_tr, lambda grid down to 1e-6.
.cv_ridge_r2_old <- function(X, y, nfold = 10L, seed = 1L) {
  n <- length(y); set.seed(seed); fold <- sample(rep_len(seq_len(nfold), n)); pred <- numeric(n)
  for (f in seq_len(nfold)) {
    tr <- fold != f; mu <- colMeans(X[tr, , drop = FALSE]); yb <- mean(y[tr])
    Xc <- sweep(X[tr, , drop = FALSE], 2, mu); yc <- y[tr] - yb; s <- svd(Xc); d2 <- s$d^2
    uy <- crossprod(s$u, yc); u2 <- s$u^2; lam <- max(d2) * 10^seq(-6, 1, length.out = 30)
    loo <- vapply(lam, function(l) { w <- d2 / (d2 + l); h <- as.numeric(u2 %*% w)
      mean(((yc - s$u %*% (w * uy)) / (1 - h))^2) }, numeric(1))
    l <- lam[which.min(loo)]; beta <- s$v %*% ((s$d / (d2 + l)) * uy)
    pred[!tr] <- yb + sweep(X[!tr, , drop = FALSE], 2, mu) %*% beta
  }
  1 - sum((y - pred)^2) / sum((y - mean(y))^2)
}
set.seed(11); nD <- 100; yn <- rnorm(nD)
for (dd in c(20, 90, 100, 200)) {
  Xn <- matrix(rnorm(nD * dd), nD); rn <- cv_ridge_r2(Xn, yn); ro <- .cv_ridge_r2_old(Xn, yn)
  cat(sprintf("noise y, n=%d d=%3d: cv R2 new %.3f (old, buggy: %.2f)\n", nD, dd, rn, ro))
  stopifnot(rn >= -0.3, rn <= 0.05)
}
X90 <- matrix(rnorm(nD * 90), nD); ys <- X90[, 1:3] %*% c(2, 2, 2) + rnorm(nD)  # beta=1 gives only ~0.28: 87 noise dims vs n_tr=90 is genuinely hard
rs <- cv_ridge_r2(X90, as.numeric(ys)); cat(sprintf("signal y, d=90: cv R2 %.3f\n", rs)); stopifnot(rs > 0.5)

## ---- CLIs ----
run(file.path(here, "03_fit_cp.R"), "combat_synth", "3", "biomni", file.path(out, "fit_biomni.rds"))
run(file.path(here, "03_fit_cp.R"), "combat_synth", "3", "repo", file.path(out, "fit_repo.rds"))
fitarg <- paste0("cp_biomni=", file.path(out, "fit_biomni.rds"), ",cp_repo=", file.path(out, "fit_repo.rds"))
run(file.path(here, "02_metrics.R"), "combat_synth", fitarg, file.path(out, "eval.csv"))
run(file.path(here, "04_corrected_eval.R"), "combat_synth", fitarg, file.path(out, "corrected.csv"))
e <- read.csv(file.path(out, "eval.csv")); cc <- read.csv(file.path(out, "corrected.csv"))
stopifnot(identical(names(e), c("representation", "covariate", "type", "metric", "value", "panel")),
          all(c("floor", "floor_masked", "metric_variant") %in% names(cc)),
          setequal(unique(e$representation), c("random", "composition", "clr_composition", "global_pca", "perct_pca", "grouped_pca", "cp_biomni", "cp_repo")),
          !anyNA(cc$floor))
## Biomni-style rows of 04 must equal 02 up to knn tie-breaking seed
b <- cc[cc$metric_variant %in% c("biomni_spearman"), ]; m <- merge(e, b, by = c("representation", "covariate"))
stopifnot(all(abs(m$value.x - m$value.y) < 1e-9))
print(unique(cc[cc$representation == "cp_biomni", c("covariate", "metric_variant", "value", "floor", "floor_masked", "n_dim")]), digits = 2)
run(file.path(here, "compare_to_biomni.R"), file.path(out, "eval.csv"))
cat("ALL SYNTHETIC CHECKS PASSED\n")
