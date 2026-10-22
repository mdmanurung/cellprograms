## smoke test on fake score/loading fits (no flashier needed for the 4 fast
## methods; the block_E_cal section at the end is guarded on flashier/ebnm):
## R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/test_methods.R
## v2: hard assertions (Defect E fix) — null must give 0 clusters for every
## method; planted must recover the R100 pattern (Defect D regression guard).
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("methods.R", "eb_cp.R", "scoring.R")) source(file.path(here, f))
set.seed(1)
cts <- sprintf("ct%02d", 1:10); ids <- sprintf("d%03d", 1:100)
strata <- setNames(sample(c("a", "b", "c"), 100, TRUE), ids)
mk <- function(planted, null = FALSE, same_genes = TRUE, k = 6, genes = 400) {
  zl <- lapply(planted, function(m) setNames(rnorm(100), ids))
  wl <- lapply(planted, function(m) rnorm(genes))
  sc <- ld <- list()
  for (ct in cts) {
    d <- ids[sort(sample(100, if (ct %in% c("ct08", "ct09", "ct10")) 60 else 100))]
    Z <- matrix(rnorm(length(d) * k), length(d), k, dimnames = list(d, paste0(ct, "_", 1:k)))
    W <- matrix(rnorm(genes * k), genes, k, dimnames = list(paste0("g", 1:genes), colnames(Z)))
    if (!null) for (p in seq_along(planted)) if (ct %in% planted[[p]]) {
      Z[, p] <- 2 * zl[[p]][d] + rnorm(length(d), sd = 0.5)
      W[, p] <- if (same_genes) wl[[p]] + rnorm(genes, sd = 0.3) else rnorm(genes)
    }
    sc[[ct]] <- Z; ld[[ct]] <- W
  }
  list(cell_types = cts, scores = sc, loadings = ld)
}
planted <- list(R100 = cts, R50 = cts[c(1, 3, 5, 8, 10)], R20 = cts[c(2, 9)])
f <- mk(planted)
show <- function(r) { cat(r$method, ":\n"); print(aggregate(cell_type ~ cluster, r$clusters, function(x) paste(sort(x), collapse = ","))) }
for (m in c("node_cor", "subspace_sum", "pa_cc", "gene_match")) show(run_method(m, f, strata, n_perm = 199, seed = 1))

chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")

## Defect D regression: node_cor must find the planted R100 pattern
r <- run_method("node_cor", f, strata, n_perm = 199, seed = 1)
big <- r$clusters$cell_type[r$clusters$cluster == which.max(table(r$clusters$cluster))]
chk(length(big) >= 8, "node_cor recovers an R100-sized cluster on planted data")

## Defect E regression: null dataset -> 0 clusters for every method
cat("--- null (expect no clusters) ---\n")
fn <- mk(planted, null = TRUE)
for (m in c("node_cor", "subspace_sum", "pa_cc", "gene_match")) {
  r <- run_method(m, fn, strata, n_perm = 199, seed = 2)
  cat(m, "clusters:", length(unique(r$clusters$cluster)), "\n")
  chk(length(unique(r$clusters$cluster)) == 0, paste(m, "returns 0 clusters on the null"))
}

## scorer self-check
score_selfcheck()

## block_E_cal: permutation-calibrated block-mass calls (needs flashier)
if (requireNamespace("flashier", quietly = TRUE) && requireNamespace("ebnm", quietly = TRUE)) {
  mk_expr_fit <- function(null = FALSE) {
    cts2 <- c("ctA", "ctB", "ctC"); ids2 <- sprintf("d%03d", 1:60)
    genes <- paste0("g", 1:200)
    z <- setNames(rnorm(60), ids2)
    fits <- list()
    for (ct in cts2) {
      Y <- matrix(rnorm(60 * 200, sd = 0.6), 60, 200, dimnames = list(ids2, genes))
      if (!null) {
        w <- rnorm(200); w[-(1:30)] <- 0
        Y <- Y + 2 * outer(z, w)          # one recurring program, all 3 types
      }
      fits[[ct]] <- list(Y_used = Y)
    }
    list(cell_types = cts2, fits = fits)
  }
  strata2 <- setNames(sample(c("a", "b"), 60, TRUE), sprintf("d%03d", 1:60))

  fpl <- mk_expr_fit(FALSE)
  fp <- run_method("block_E_cal", fpl, strata2, n_perm = 199, seed = 3)
  chk(nrow(fp$clusters) >= 1L, "block_E_cal keeps the planted call")
  chk(all(fp$clusters$fwer_p <= 0.05), "block_E_cal fills significant fwer_p on planted calls")
  chk(setequal(fp$clusters$cell_type[fp$clusters$cluster == fp$clusters$cluster[1]],
               c("ctA", "ctB", "ctC")), "block_E_cal planted call spans all 3 types")
  chk(all(fp$clusters$stat > 0 & fp$clusters$stat <= 1),
      "block_E_cal stat is the mean pairwise |cor|")

  fnn <- mk_expr_fit(TRUE)
  be <- run_method("block_E", fnn, strata2, seed = 4)
  cat("block_E on null calls:", length(unique(be$clusters$cluster)), "clusters\n")
  rn <- run_method("block_E_cal", fnn, strata2, n_perm = 199, seed = 4, bfit = be$fit)
  chk(nrow(rn$clusters) == 0L, "block_E_cal returns 0 clusters on the null")

  fp2 <- run_method("block_E_cal", fpl, strata2, n_perm = 199, seed = 3, bfit = fp$fit)
  chk(identical(fp$clusters, fp2$clusters), "block_E_cal cached-bfit run reproduces clusters")
}

cat("all method smoke checks passed\n")
