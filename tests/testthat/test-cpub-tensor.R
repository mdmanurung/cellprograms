test_that("CPUB-001 tensor normalization preserves products for either anchor sign", {
  for (sgn in c(-1, 1)) {
    u <- sgn * seq_len(10)
    v <- c(-2, 0, 3, 1)
    w <- c(1, -2, 0)
    n <- .cp_canonical_component(u, v, w)
    expect_gte(n$u[which.max(abs(n$u))], 0)
    for (ci in seq_along(w)) {
      expect_equal(n$w[ci] * tcrossprod(n$v, n$u),
                   w[ci] * tcrossprod(v, u), tolerance = 1e-10)
    }
    expect_equal(sqrt(mean(n$v^2)), 1)
    expect_identical(n$w, w)
  }
  n <- .cp_canonical_component(rep(0, 10), rep(0, 4), c(1, 2))
  expect_true(all(is.finite(unlist(n))))
  # Independent two-component sum oracle, not just one normalized factor.
  set.seed(13)
  U <- matrix(rnorm(30), 10, 3); V <- matrix(rnorm(12), 4, 3)
  W <- matrix(rnorm(6), 2, 3)
  normal <- lapply(1:3, function(k) .cp_canonical_component(U[, k], V[, k], W[, k]))
  for (ct in 1:2) {
    before <- after <- matrix(0, 4, 10)
    for (k in 1:3) {
      before <- before + W[ct, k] * tcrossprod(V[, k], U[, k])
      after <- after + normal[[k]]$w[ct] * tcrossprod(normal[[k]]$v, normal[[k]]$u)
    }
    expect_equal(before, after, tolerance = 1e-12)
  }
})

test_that("CPUB-001 actual assembly and diagnostic agree with mocked rank-one solver", {
  ids <- paste0("d", 1:20)
  genes <- paste0("g", 1:10)
  u <- -seq_len(10)/10
  v <- setNames(seq_len(20)/10, ids)
  w <- c(1, -2)
  mats <- lapply(w, function(wi) {
    m <- wi * outer(v, u); dimnames(m) <- list(ids, genes); m
  })
  names(mats) <- c("A", "B")
  mats$B <- mats$B[-c(2, 4), , drop = FALSE]
  scope <- new.env(parent = environment(fit_eb_cp))
  scope$.cp_rank1 <- function(R, obs, tau2, ...) {
    list(u = u, v = unname(v[obs$donors]), w = w, tau2 = tau2)
  }
  scope$.cp_heldout_z <- function(...) 100
  local_fit <- fit_eb_cp
  environment(local_fit) <- scope
  cp <- local_fit(mats, Kmax = 1, backfit_iter = 0, seed = 3)
  ss <- total <- 0
  for (ct in names(mats)) {
    pred <- cp$W[ct, 1] * tcrossprod(cp$V[rownames(mats[[ct]]), 1], cp$U[, 1])
    expect_equal(unname(pred), unname(mats[[ct]]), tolerance = 1e-10)
    ss <- ss + sum((pred - mats[[ct]])^2)
    total <- total + sum(mats[[ct]]^2)
  }
  expect_equal(cp$frac_unexplained, ss/total, tolerance = 1e-10)
})

test_that("CPUB-001 zero-rank output is well-formed and zero-energy error undefined", {
  m <- matrix(0, 20, 10, dimnames = list(paste0("d", 1:20), paste0("g", 1:10)))
  cp <- fit_eb_cp(list(A = m, B = m), Kmax = 0)
  expect_equal(dim(cp$U), c(10, 0))
  expect_equal(dim(cp$V), c(20, 0))
  expect_equal(dim(cp$W), c(2, 0))
  expect_true(is.na(cp$frac_unexplained))
})

test_that("CPUB-009 unusable strata and ragged halves fail before solving", {
  ids <- paste0("d", 1:30)
  expect_error(.cp_split_donors(ids, setNames(1:30, ids), 1), "singleton/continuous")
  expect_error(.cp_split_donors(ids, rep("a", 30), 1), "unique donor names")
  m <- matrix(seq_len(300), 30, 10, dimnames = list(ids, paste0("g", 1:10)))
  expect_error(fit_eb_cp(list(A = m, B = m[1:2, , drop = FALSE]), Kmax = 0),
               "both EB-CP split halves")
})

test_that("CPUB-001 tiny real EB normal-means fit has self-consistent diagnostics", {
  skip_if_not_installed("ebnm")
  set.seed(7)
  ids <- paste0("d", 1:30); genes <- paste0("g", 1:20)
  z <- rnorm(30); loading <- c(rep(-1, 10), rep(1, 10))
  mats <- lapply(c(1, -1), function(w) {
    y <- w * outer(z, loading) + matrix(rnorm(600, sd = .1), 30)
    dimnames(y) <- list(ids, genes); y
  })
  names(mats) <- c("A", "B")
  cp <- fit_eb_cp(mats, Kmax = 2, heldout_z = -1e6, max_iter = 10,
                  backfit_iter = 1, seed = 3)
  err <- total <- 0
  for (ct in names(mats)) {
    pred <- matrix(0, 30, 20)
    for (k in seq_len(cp$n_components)) {
      pred <- pred + cp$W[ct, k] * tcrossprod(cp$V[ids, k], cp$U[, k])
    }
    err <- err + sum((mats[[ct]] - pred)^2)
    total <- total + sum(mats[[ct]]^2)
  }
  expect_equal(cp$frac_unexplained, err/total, tolerance = 1e-10)
})
