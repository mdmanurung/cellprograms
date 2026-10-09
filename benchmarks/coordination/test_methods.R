## smoke test on fake score/loading fits (no flashier): R_LIBS_USER=/nonexistent Rscript benchmarks/coordination/test_methods.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats); library(utils)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
source(file.path(here, "methods.R"))
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
show <- function(r) {
  cat(r$method, ":\n")
  if (nrow(r$clusters)) print(aggregate(cell_type ~ cluster, r$clusters, function(x) paste(sort(x), collapse = ",")))
  else cat("  (no clusters)\n")
}
M4 <- c("node_cor", "subspace_sum", "pa_cc", "gene_match")
for (m in M4) show(run_method(m, f, strata, n_perm = 199, seed = 1))

## .perm_rows must move values, keep row names and stay inside strata (a no-op permutation made node_cor never reject)
Mx <- f$scores$ct01; st1 <- strata[rownames(Mx)]
set.seed(5); Mp <- .perm_rows(Mx, st1)
stopifnot(identical(rownames(Mp), rownames(Mx)), !isTRUE(all.equal(unname(Mp), unname(Mx))))
for (s in unique(st1)) stopifnot(identical(sort(unname(Mp[st1 == s, 1])), sort(unname(Mx[st1 == s, 1]))))

## null: every method must return 0 clusters (regression guard for no-op permutations)
cat("--- null (expect no clusters) ---\n")
fn <- mk(planted, null = TRUE)
for (m in M4) {
  r <- run_method(m, fn, strata, n_perm = 199, seed = 2)
  n <- length(unique(r$clusters$cluster)); cat(m, "clusters:", n, "\n")
  stopifnot(n == 0L)
}
## power: with 3 planted programs node_cor must recover at least one multi-cell-type cluster
rp <- run_method("node_cor", f, strata, n_perm = 199, seed = 1)
stopifnot(nrow(rp$clusters) > 0L)
cat("OK\n")
