.fake_fit <- function(Za, Zb) {
  structure(list(scores = list(A = Za, B = Zb), loadings = list(A = Za[0, ], B = Zb[0, ])),
            class = "cell_program_fit")
}
.rand_scores <- function(n, k, ct) {
  matrix(rnorm(n * k), n, k, dimnames = list(paste0("s", seq_len(n)), paste0(ct, "_", seq_len(k))))
}

test_that("pair-level permutation p is calibrated under independence and powered under sharing", {
  set.seed(42)
  n <- 40; reps <- 150
  p0 <- p1 <- numeric(reps)
  for (r in seq_len(reps)) {
    Za <- .rand_scores(n, 3, "A"); Zb <- .rand_scores(n, 3, "B")
    p0[r] <- principal_angles(.fake_fit(Za, Zb), "A", "B", n_perm = 49, seed = r)$p_pair
    Zb[, 1] <- Za[, 1] * 0.8 + rnorm(n, sd = 0.6)   # one shared program
    p1[r] <- principal_angles(.fake_fit(Za, Zb), "A", "B", n_perm = 49, seed = r)$p_pair
  }
  expect_lt(mean(p0 < 0.05), 0.12)          # nominal 0.05; binomial slack at 150 reps
  expect_gt(mean(p0 < 0.05), 0.005)         # not degenerate-conservative
  expect_gt(mean(p1 < 0.05), 0.9)
})

test_that("sharing_spectrum applies BH over callable pairs and flags underpowered ones", {
  set.seed(1)
  n <- 30
  mk <- function(k, ct) .rand_scores(n, k, ct)
  Z1 <- mk(3, "A"); Z2 <- mk(3, "B"); Z3 <- mk(2, "C")
  Z2[, 1] <- Z1[, 1] + rnorm(n, sd = 0.3)             # A-B truly shared
  Z4 <- mk(20, "D"); Z5 <- mk(20, "E")                # sqrt(20)*2/sqrt(30) > 1: uncallable
  fit <- structure(list(scores = list(A = Z1, B = Z2, C = Z3, D = Z4, E = Z5),
                        loadings = list(A = Z1[0, ], B = Z2[0, ], C = Z3[0, ], D = Z4[0, ], E = Z5[0, ])),
                   class = "cell_program_fit")
  ss <- sharing_spectrum(fit, pairs = list(c("A", "B"), c("A", "C"), c("D", "E")), n_perm = 99, seed = 1)$summary
  expect_equal(ss$underpowered, c(FALSE, FALSE, TRUE))
  expect_true(is.na(ss$q_pair[3]) && !ss$shared[3])
  expect_true(ss$shared[1])
  expect_false(ss$shared[2])
  expect_equal(ss$q_pair[1:2], stats::p.adjust(ss$p_pair[1:2], "BH"))
})
