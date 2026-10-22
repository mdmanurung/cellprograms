here <- "benchmarks/coordination"; root <- normalizePath(".")
library(stats)
for (f in list.files("R", full.names = TRUE)) source(f)
for (f in c("sim_coord.R", "eb_cp.R")) source(file.path(here, f))
options(ms.G = 1200L)

## instrumented rank-1: report gene-mode and ct-mode log-marginal-likelihood gains
probe <- function(mode, seed, amp = 0.75, breadth = "R100") {
  d <- gen_coord(mode, amp = amp, seed = seed, breadth = breadth)
  mats <- adjusted_matrices(d$x, "institute")
  genes <- Reduce(intersect, lapply(mats, colnames))
  mats <- lapply(mats, function(m) m[, genes, drop = FALSE])
  donors <- sort(unique(unlist(lapply(mats, rownames))))
  obs <- list(donors = donors, idx = lapply(mats, function(m) match(rownames(m), donors)))
  tau2 <- vapply(mats, function(m) max(mean(m^2), 1e-8), numeric(1))
  ## one rank-1 step, tracking LL gains at the final iteration
  R <- mats
  f <- local({
    C <- length(R); G <- ncol(R[[1]]); D <- length(donors)
    init <- .cp_seed(R, obs); v <- init$v; w <- init$w; u <- init$u
    gains <- NULL
    for (it in 1:60) {
      num <- numeric(G); den <- 0; s2 <- 0
      for (ci in seq_len(C)) {
        vv <- v[obs$idx[[ci]]]
        num <- num + w[ci] * drop(crossprod(vv, R[[ci]]))
        s2v <- sum(vv^2); den <- den + w[ci]^2 * s2v; s2 <- s2 + w[ci]^2 * s2v * tau2[ci]
      }
      x <- num / den; s <- sqrt(s2) / den
      eg <- ebnm::ebnm_point_laplace(x, s)
      ll_gene <- eg$log_likelihood - sum(dnorm(x, 0, s, log = TRUE))
      u <- as.numeric(eg$posterior$mean); u2 <- sum(u^2)
      num <- numeric(D); den <- numeric(D); s2 <- numeric(D)
      for (ci in seq_len(C)) {
        ii <- obs$idx[[ci]]; proj <- drop(R[[ci]] %*% u)
        num[ii] <- num[ii] + w[ci] * proj; den[ii] <- den[ii] + w[ci]^2 * u2
        s2[ii] <- s2[ii] + w[ci]^2 * u2 * tau2[ci]
      }
      ok <- den > 1e-12
      ed <- ebnm::ebnm_normal(num[ok] / den[ok], sqrt(s2[ok]) / den[ok])
      v[ok] <- as.numeric(ed$posterior$mean); v[!ok] <- 0
      num <- numeric(C); den <- numeric(C); s2 <- numeric(C)
      for (ci in seq_len(C)) {
        vv <- v[obs$idx[[ci]]]; s2v <- sum(vv^2)
        num[ci] <- sum(vv * drop(R[[ci]] %*% u)); den[ci] <- s2v * u2; s2[ci] <- tau2[ci] * s2v * u2
      }
      okc <- den > 1e-12
      ec <- ebnm::ebnm_point_laplace(num[okc] / den[okc], sqrt(s2[okc]) / den[okc])
      ll_ct <- ec$log_likelihood - sum(dnorm(num[okc] / den[okc], 0, sqrt(s2[okc]) / den[okc], log = TRUE))
      w[okc] <- as.numeric(ec$posterior$mean); w[!okc] <- 0
      ss <- 0
      for (ci in seq_len(C)) {
        E <- R[[ci]] - w[ci] * tcrossprod(v[obs$idx[[ci]]], u)
        tau2[ci] <- max(mean(E^2), 1e-8); ss <- ss + sum(E^2)
      }
      gains <- c(ll_gene = ll_gene, ll_ct = ll_ct)
      if (it > 1 && abs(gains_prev - ss) < 1e-6 * (1 + abs(ss))) break
      gains_prev <- ss
    }
    list(u = u, v = v, w = w, gains = gains, it = it)
  })
  data.frame(mode = mode, seed = seed, ll_gene = f$gains["ll_gene"], ll_ct = f$gains["ll_ct"],
             max_abs_w = max(abs(f$w)), rms_u = sqrt(mean(f$u^2)), it = f$it)
}

res <- rbind(
  do.call(rbind, lapply(1:8, function(s) probe("null", 200 + s))),
  do.call(rbind, lapply(1:4, function(s) probe("act+genes", 300 + s, breadth = "R100"))),
  do.call(rbind, lapply(1:4, function(s) probe("act_only", 400 + s, amp = 0.5, breadth = "R50")))
)
print(res, digits = 3)
