# Workstream E extension: detection scenario under dropout/zero-inflation.
# Underlying data identical to the clean detection datasets (same seed stream
# in gen_detection), with entrywise Bernoulli dropout applied to the COUNT
# matrices. Masks are drawn once per rep per cell type (seeds 6000+rep,
# 7000+rep) and are nested across rates (dropout iff u < rate) and shared
# across rho (pairing across rho preserved).
n <- 80L; p <- 600L

gen_detection <- function(seed, rho) {
  set.seed(1000 + seed)
  k <- 10L
  ids <- paste0("d", seq_len(n))
  site <- rep(c("s1", "s2"), length.out = n)
  z <- rnorm(n)
  scores <- list()
  for (ct in c("A", "B")) {
    S <- matrix(rnorm(n * k), n, k, dimnames = list(ids, paste0(ct, "_", seq_len(k))))
    if (rho > 0) S[, 1] <- sqrt(rho) * z + sqrt(1 - rho) * rnorm(n)
    scores[[ct]] <- S
  }
  Y <- lapply(scores, function(S) {
    M <- S %*% t(matrix(rnorm(p * k, sd = 0.1), p, k)) +
      matrix(rnorm(n * p, 5, 0.5), n, p, dimnames = list(ids, paste0("G", seq_len(p))))
    dimnames(M) <- list(ids, paste0("G", seq_len(p)))
    M
  })
  list(Y = Y, z = z, rho = rho, site = site, scores = scores)
}

for (rep in 1:20) {
  set.seed(6000 + rep); u_A <- runif(n * p)
  set.seed(7000 + rep); u_B <- runif(n * p)
  for (rate in c(0.1, 0.3)) {
    for (rho in c(0, 0.3, 0.45, 0.6)) {
      d <- gen_detection(rep, rho)
      tag <- sprintf("det_drop%g_rho%g_%02d", rate, rho, rep)
      for (ct in c("A", "B")) {
        counts <- round(exp(d$Y[[ct]]))
        u <- if (ct == "A") u_A else u_B
        counts[u < rate] <- 0
        write.csv(counts, sprintf("/workspace/audit/cnmf_bench/data/%s_%s_counts.csv", tag, ct), row.names = TRUE)
      }
      saveRDS(list(z = d$z, rho = rho, site = d$site, rate = rate, scores = d$scores),
              sprintf("/workspace/audit/cnmf_bench/data/%s.rds", tag))
    }
  }
}
cat("detection dropout scenarios generated\n")
