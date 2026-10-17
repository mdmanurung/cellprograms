# Workstream E extension: dropout/zero-inflation recovery scenarios.
# Same latent data as recovery_k3 (identical seed stream), with entrywise
# Bernoulli dropout applied to the COUNT matrix. Masks are nested across
# rates (u drawn once per rep; dropout iff u < rate). Truth (z, W) unchanged.
# Note: uniform entrywise dropout is an adversarial-to-both-methods stress
# test; real pseudobulk zeros are mean-dependent, which is not modeled here.
dir.create("/workspace/audit/cnmf_bench/data", showWarnings = FALSE)
n <- 80L; p <- 600L; block <- 30L

gen_recovery <- function(seed, k_true) {
  set.seed(seed)
  z <- matrix(rnorm(n * k_true), n, k_true)
  W <- matrix(0, p, k_true, dimnames = list(paste0("G", seq_len(p)), paste0("P", seq_len(k_true))))
  for (k in seq_len(k_true)) W[((k - 1) * block + 1):(k * block), k] <- 1
  Y <- matrix(rnorm(n * p, 5, 0.5), n, p, dimnames = list(paste0("d", seq_len(n)), paste0("G", seq_len(p))))
  Y <- Y + (z + rnorm(n * k_true, sd = 0.3)) %*% t(W)
  list(Y = Y, z = z, W = W, k_true = k_true)
}

for (rep in 1:10) {
  d <- gen_recovery(rep, 3L)
  set.seed(5000 + rep)
  u <- runif(n * p)
  for (rate in c(0.1, 0.3)) {
    counts <- round(exp(d$Y))
    counts[u < rate] <- 0
    tag <- sprintf("recovery_drop%g_%02d", rate, rep)
    write.csv(counts, sprintf("/workspace/audit/cnmf_bench/data/%s_counts.csv", tag), row.names = TRUE)
    write.csv(d$z, sprintf("/workspace/audit/cnmf_bench/data/%s_scores.csv", tag), row.names = TRUE)
    saveRDS(list(Y = d$Y, z = d$z, W = d$W, k_true = 3L, rate = rate,
                 drop_mask = (u < rate)),
            sprintf("/workspace/audit/cnmf_bench/data/%s.rds", tag))
  }
}
cat("dropout scenarios generated\n")
