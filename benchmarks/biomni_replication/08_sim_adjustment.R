## 08_sim_adjustment.R <rep_from> <rep_to> <out_csv>
## Simulation harness: does covariate adjustment keep the sharing test calibrated (type-I under no sharing)
## while still recovering a true shared program, when a batch shifts both cell types?
## Scenarios: shared in {FALSE, TRUE} x rho in {0, 0.5} (rho = corr of the shared/own program with batch).
## Arms: none | hard (covariates=, residualize) | soft (covariates=, covariate_mode = "fixed", SOFA-style).
a <- commandArgs(TRUE); stopifnot(length(a) == 3)
reps <- seq(as.integer(a[1]), as.integer(a[2])); out <- a[3]
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
library(stats); library(utils)
for (f in list.files(file.path(here, "..", "..", "R"), full.names = TRUE)) source(f)

n <- 100; g <- 300; K <- 5; n_perm <- 499
arms <- list(none = list(NULL, "residualize"), hard = list("batch", "residualize"),
             soft = list("batch", "fixed"))

sim_one <- function(shared, rho) {
  ids <- paste0("d", seq_len(n))
  batch <- rbinom(n, 1, 0.15)                          # minor site, like St_Georges
  bz <- as.numeric(scale(batch))
  mk_prog <- function() rho * bz + sqrt(1 - rho^2) * rnorm(n)
  z <- mk_prog()                                       # the shared program (if shared)
  load <- function() { w <- numeric(g); w[sample(g, 30)] <- rnorm(30, sd = 1.5); w }
  bload <- rnorm(g, sd = 1) * (runif(g) < 0.5)         # batch hits the same genes in both cell types
  mk_ct <- function() {
    zc <- if (shared) z else mk_prog()
    Y <- outer(zc, load()) + outer(bz, bload) * 1.5 + matrix(rnorm(n * g), n, g)
    dimnames(Y) <- list(ids, paste0("g", seq_len(g)))
    list(Y = Y, z = zc)
  }
  A <- mk_ct(); B <- mk_ct()
  list(x = as_cell_program_data(list(A = A$Y, B = B$Y),
                                sample_metadata = data.frame(observation_id = ids, batch = batch)),
       z = z, zA = A$z, zB = B$z, batch = batch)
}

best_cor <- function(S, v) if (ncol(S)) max(abs(cor(S, v))) else NA_real_

rows <- list()
for (rep in reps) for (shared in c(FALSE, TRUE)) for (rho in c(0, 0.5)) {
  set.seed(rep * 100 + shared * 10 + rho * 4)
  d <- sim_one(shared, rho)
  for (arm in names(arms)) {
    fit <- canonicalize_programs(fit_celltype_programs(d$x, max_factors = K, covariates = arms[[arm]][[1]],
                                                       covariate_mode = arms[[arm]][[2]],
                                                       features = "all", seed = 1))
    SA <- fit$scores$A; SB <- fit$scores$B
    pa <- if (ncol(SA) && ncol(SB)) principal_angles(fit, "A", "B", n_perm = n_perm, seed = 1) else NULL
    rows[[length(rows) + 1]] <- data.frame(
      rep = rep, shared = shared, rho = rho, arm = arm,
      k_a = ncol(SA), k_b = ncol(SB),
      p_pair = if (is.null(pa)) NA_real_ else pa$p_pair,
      underpowered = if (is.null(pa)) NA else pa$underpowered,
      max_cor_batch = max(best_cor(SA, d$batch), best_cor(SB, d$batch), na.rm = TRUE),
      best_cor_z = mean(c(best_cor(SA, d$zA), best_cor(SB, d$zB)), na.rm = TRUE))
  }
  write.csv(do.call(rbind, rows), out, row.names = FALSE)
}
