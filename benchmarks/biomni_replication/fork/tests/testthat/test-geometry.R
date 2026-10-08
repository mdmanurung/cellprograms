## Geometry: principal angles, rotation invariance, union/intersection equivalence

test_that("principal angles are invariant to within-subspace rotation", {
  f <- get_sim_fit(4L)$fit
  pa1 <- principal_angles(f, "CD8", "NK", space = "scores", n_perm = 0L)
  Z <- f$scores$CD8
  Q <- qr.Q(qr(matrix(stats::rnorm(ncol(Z)^2), ncol(Z))))
  f2 <- f
  f2$scores$CD8 <- Z %*% Q   # same column subspace, different basis
  pa2 <- principal_angles(f2, "CD8", "NK", space = "scores", n_perm = 0L)
  expect_equal(pa1$cosines, pa2$cosines, tolerance = 1e-6)
})

test_that("disjoint subspaces give 90 deg; identical subspaces give 0 deg", {
  set.seed(1)
  N <- 80L
  A <- matrix(stats::rnorm(N * 2), N, 2)
  B_same <- A %*% matrix(c(1, 1, -1, 1), 2, 2)      # same subspace
  B_orth <- matrix(stats::rnorm(N * 2), N, 2)
  B_orth <- qr.Q(qr(cbind(A, B_orth)))[, 3:4]        # orthogonal to A
  rownames(A) <- rownames(B_same) <- rownames(B_orth) <- paste0("o", 1:N)
  fake <- list(scores = list(A = A, B = B_same, C = B_orth),
               loadings = list(A = matrix(1, 5, 2), B = matrix(1, 5, 2),
                               C = matrix(1, 5, 2)))
  pa_same <- principal_angles(fake, "A", "B", space = "scores", n_perm = 0L)
  pa_orth <- principal_angles(fake, "A", "C", space = "scores", n_perm = 0L)
  expect_equal(pa_same$angles, c(0, 0), tolerance = 1e-6)
  expect_equal(pa_orth$angles, c(pi / 2, pi / 2), tolerance = 1e-6)
})

test_that("loadings-space angles use the gene intersection (private genes ignored)", {
  set.seed(2)
  g_a <- paste0("g", 1:100); g_b <- paste0("g", 51:150)
  Wa <- matrix(stats::rnorm(200), 100, 2, dimnames = list(g_a, NULL))
  Wb <- matrix(stats::rnorm(200), 100, 2, dimnames = list(g_b, NULL))
  fake <- list(loadings = list(A = Wa, B = Wb),
               scores = list(A = matrix(1, 30, 2), B = matrix(1, 30, 2)))
  pa_pkg <- suppressWarnings(principal_angles(fake, "A", "B", space = "loadings",
                                              n_perm = 0L))
  gi <- intersect(g_a, g_b)
  pa_inter <- suppressWarnings(principal_angles(
    list(loadings = list(A = Wa[gi, ], B = Wb[gi, ]),
         scores = fake$scores),
    "A", "B", space = "loadings", n_perm = 0L))
  expect_equal(pa_pkg$cosines, pa_inter$cosines, tolerance = 1e-8)
  expect_equal(unname(pa_pkg$dims["effective_dim"]), length(gi))
})

test_that("canonicalization preserves reconstruction and unit-RMS scores", {
  sim <- simulate_programs(scenario = 4L, N = 60L, G = 150L, seed = 3L)
  f0 <- fit_celltype_programs(sim$x, max_factors = 8L)
  f1 <- canonicalize_programs(f0)
  for (ct in names(f0$scores)) {
    if (ncol(f0$scores[[ct]]) == 0L) next
    R0 <- f0$scores[[ct]] %*% t(f0$loadings[[ct]])
    R1 <- f1$scores[[ct]] %*% t(f1$loadings[[ct]])
    expect_equal(R0, R1, tolerance = 1e-6)
    rms <- sqrt(colMeans(f1$scores[[ct]]^2))
    expect_equal(unname(rms), rep(1, ncol(f1$scores[[ct]])), tolerance = 1e-6)
  }
})
