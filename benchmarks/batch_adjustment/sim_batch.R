# batch-adj-v1 simulator extension (see PREREG.md in this directory).
# Plants gene-level batch location/scale shifts on simulate_pseudobulk output,
# and realizes count-level data for the ComBat-seq arm.

# Plant batch effects: B batches (balanced, shuffled), a shared shifted gene
# subset per batch (same chemistry across cell types), plus a per-batch global
# log-scale shift. `confound_a` != 0 makes batch quantiles track a planted
# program's donor activities (secondary confounded arm; a = target correlation).
plant_batch <- function(sim, n_batches = 4, sigma_b = 0.375,
                        shift_fraction = 0.10, scale_sd = 0.05,
                        confound_a = 0, seed = 2) {
  set.seed(seed)
  pb <- sim$pseudobulk
  n <- nrow(pb[[1]])
  genes <- colnames(pb[[1]])
  if (confound_a > 0) {
    ct0 <- names(sim$truth)[1]
    z0 <- scale(sim$truth[[ct0]][[1]]$z)[, 1]
    noise <- rnorm(n, 0, 1)
    score <- confound_a * z0 + sqrt(1 - confound_a^2) * noise
    brk <- seq(0, 1, length.out = n_batches + 1)
    batch <- factor(paste0("b", cut(rank(score), breaks = brk * n,
                                    labels = seq_len(n_batches), include.lowest = TRUE)))
  } else {
    batch <- factor(sample(paste0("b", rep(seq_len(n_batches), length.out = n))))
  }
  shift_genes <- sample(length(genes), max(1L, round(shift_fraction * length(genes))))
  for (b in seq_len(n_batches)) {
    idx <- which(batch == paste0("b", b))
    delta <- rep(0, length(genes))
    delta[shift_genes] <- rnorm(length(shift_genes), 0, sigma_b)
    scale_shift <- rnorm(1, 0, scale_sd)
    for (ct in names(pb)) {
      M <- pb[[ct]][idx, , drop = FALSE]
      pb[[ct]][idx, ] <- sweep(M, 2, delta, "+") + scale_shift
    }
  }
  list(pseudobulk = pb,
       batch = setNames(batch, rownames(pb[[1]])),
       shift_genes = genes[shift_genes],
       sigma_b = sigma_b, confound_a = confound_a)
}

# Count-level realization for the ComBat-seq arm: treat simulated values as
# log-expression, exponentiate with donor library-size factors, round.
make_counts <- function(pb, seed = 3) {
  set.seed(seed)
  lapply(pb, function(Y) {
    lib <- exp(rnorm(nrow(Y), 0, 0.4))
    M <- sweep(exp(Y), 1, lib, "*")
    M <- matrix(round(M), nrow(M), ncol(M), dimnames = dimnames(Y))
    M
  })
}

# Standard log1p-CPM normalization of a counts matrix (donors x genes).
lognorm_counts <- function(counts) {
  lib <- colSums(counts)
  t(log1p(t(sweep(counts, 2, lib, "/")) * 1e4))
}
