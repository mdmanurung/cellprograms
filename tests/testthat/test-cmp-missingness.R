# Regression tests for missing-donor-score handling in CMP detection (audit finding B-6).
#
# Pre-fix behavior (repo @ 413e290): rank()'s default na.last=TRUE converted NA
# scores into top ranks, correlations were computed over the corrupted scores,
# and edge `n` reported the nominal donor count. A planted pair with 20%
# missingness was silently missed. Post-fix contract:
#   - NA scores stay NA through rank transformation and residualization;
#   - cross-cell-type correlations use pairwise-complete donors;
#   - edge `n` reports the pairwise-complete donor count for that program pair;
#   - find_coordinated_programs() warns when any tested program has missing
#     donor scores, and is silent on complete data.

# Toy canonicalized fit (same generator as test-coordination.R, duplicated here
# because test files run alphabetically and cannot rely on each other's objects).
.cmp_toy <- function(n = 80, k = 3, g = 200, seed = 1, couple = TRUE, site_shift = 0, context_only = FALSE,
                     cts = c("A", "B", "C")) {
  set.seed(seed)
  ids <- sprintf("d%03d", seq_len(n))
  site <- rep(c("s1", "s2"), length.out = n)
  ctx <- rep(c(0, 1), each = n / 2)
  z <- rnorm(n)
  genes <- sprintf("G%03d", seq_len(g))
  scores <- loadings <- pb <- list()
  for (ct in cts) {
    S <- matrix(rnorm(n * k), n, k, dimnames = list(ids, paste0(ct, "_", seq_len(k))))
    S <- S + site_shift * (site == "s1")
    if (couple && ct %in% cts[1:3]) {
      j <- match(ct, cts)
      S[, j] <- (if (context_only) 2 * ctx else z) * (if (ct == cts[3]) -1 else 1) + rnorm(n, sd = 0.6)
    }
    W <- matrix(rnorm(g * k, sd = 0.1), g, k, dimnames = list(genes, colnames(S)))
    W[(1:20) + 20 * (match(ct, cts) - 1), ] <- W[(1:20) + 20 * (match(ct, cts) - 1), ] + 1
    scores[[ct]] <- S; loadings[[ct]] <- W
    pb[[ct]] <- matrix(rnorm(n * g, 5), n, g, dimnames = list(ids, genes))
  }
  structure(list(scores = scores, loadings = loadings, cell_types = cts, observation_ids = ids, pseudobulk = pb,
                 sample_metadata = data.frame(observation_id = ids, site = site, ctx = ctx, stringsAsFactors = FALSE)),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

# inject NA donors into every program of the given cell types (per-type pattern)
.na_by_type <- function(fit, types, n_miss = 16L, seed = 42) {
  set.seed(seed)
  for (ct in types) {
    S <- fit$scores[[ct]]
    S[sample(nrow(S), n_miss), ] <- NA_real_
    fit$scores[[ct]] <- S
  }
  fit
}

# inject NA donors into each program independently (per-program pattern)
.na_by_program <- function(fit, n_miss = 16L, seed = 42) {
  set.seed(seed)
  for (ct in names(fit$scores)) {
    S <- fit$scores[[ct]]
    for (j in seq_len(ncol(S))) S[sample(nrow(S), n_miss), j] <- NA_real_
    fit$scores[[ct]] <- S
  }
  fit
}

test_that("rank-gauss keeps NA scores as NA instead of top ranks", {
  S <- matrix(c(1, 3, 2, NA, NA, 1, 4, 5, 5, 2), 5, 2, byrow = TRUE,
              dimnames = list(paste0("d", 1:5), c("p1", "p2")))
  R <- cellprograms:::.cmp_rank_gauss(S)
  expect_true(is.na(R["d3", "p1"]))
  expect_true(is.na(R["d2", "p2"]))
  expect_false(anyNA(R[-3, "p1"]))
  expect_false(anyNA(R[-2, "p2"]))
  # observed values are proper normal scores: rank 1 of 4 observed -> qnorm(0.5/4)
  expect_equal(R["d1", "p1"], stats::qnorm(0.5 / 4))
  expect_equal(R["d5", "p1"], stats::qnorm(3.5 / 4))
})

test_that("residualization keeps NA rows NA and residualizes observed rows", {
  S <- matrix(rnorm(50), 10, 5, dimnames = list(paste0("d", 1:10), paste0("p", 1:5)))
  S[1:2, 1] <- NA_real_
  Z <- cellprograms:::.cmp_prep(list(A = S), setNames(rep("all", 10), rownames(S)), min_donors = 5)
  expect_true(all(is.na(Z$A[1:2, 1])))
  expect_false(anyNA(Z$A[-(1:2), 1]))
  expect_false(anyNA(Z$A[, -1]))
  # observed rows are residualized: mean ~ 0
  expect_lt(abs(mean(Z$A[-(1:2), 1])), 1e-10)
})

test_that("edge statistics use pairwise-complete donors and honest n", {
  n <- 80
  ids <- sprintf("d%03d", seq_len(n))
  Za <- matrix(rnorm(n * 2), n, 2, dimnames = list(ids, c("A_1", "A_2")))
  Zb <- matrix(rnorm(n * 2), n, 2, dimnames = list(ids, c("B_1", "B_2")))
  Za[1:16, "A_1"] <- NA_real_            # per-program: only A_1 has missing donors
  Zb[9:24, "B_1"] <- NA_real_            # overlapping donors 9-16 also missing in B_1
  Z <- list(A = Za, B = Zb)
  pairs <- list(list(a = "A", b = "B", ids = ids))
  E <- cellprograms:::.cmp_edges(Z, pairs)
  n_ref <- function(pa, pb) length(intersect(
    ids[is.finite(Za[, pa])], ids[is.finite(Zb[, pb])]))
  expect_equal(E$n[E$prog_a == "A_1" & E$prog_b == "B_1"], n_ref("A_1", "B_1"))  # 56
  expect_equal(E$n[E$prog_a == "A_1" & E$prog_b == "B_2"], n_ref("A_1", "B_2"))  # 64
  expect_equal(E$n[E$prog_a == "A_2" & E$prog_b == "B_1"], n_ref("A_2", "B_1"))  # 64
  expect_equal(E$n[E$prog_a == "A_2" & E$prog_b == "B_2"], n_ref("A_2", "B_2"))  # 80
  # correlation matches a manual pairwise-complete computation
  keep <- intersect(ids[is.finite(Za[, "A_1"])], ids[is.finite(Zb[, "B_1"])])
  r_manual <- cor(Za[keep, "A_1"], Zb[keep, "B_1"])
  e <- E[E$prog_a == "A_1" & E$prog_b == "B_1", ]
  expect_equal(e$stat, abs(atanh(r_manual) * sqrt(length(keep) - 3)), tolerance = 1e-10)
  # complete data: n equals the nominal donor count (backward compatible)
  Zc <- lapply(Z, function(M) { M[] <- rnorm(length(M)); M })
  E2 <- cellprograms:::.cmp_edges(Zc, pairs)
  expect_true(all(E2$n == n))
})

test_that("find_coordinated_programs warns on missing scores and is silent on complete data", {
  fit <- .cmp_toy()
  expect_no_warning(res0 <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none"))
  fit_m <- .na_by_type(.cmp_toy(), "A")
  expect_warning(res <- find_coordinated_programs(fit_m, strata = "site", n_perm = 99, ambient = "none"),
                 regexp = "pairwise-complete")
})

test_that("planted pair is detected under per-type missingness (20% in one cell type)", {
  fit <- .na_by_type(.cmp_toy(), "A", n_miss = 16L, seed = 11)
  # warning contract is verified in its own block; suppress here for a clean run
  res <- suppressWarnings(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none", seed = 11))
  expect_true(all(res$cmps$node_p < 0.05))
  expect_setequal(paste(res$cmps$cell_type, res$cmps$program_id), c("A A_1", "B B_2", "C C_3"))
  # honest edge n: A-B edges use 64 pairwise-complete donors, not 80
  ab <- res$edges[res$edges$ct_a == "A" | res$edges$ct_b == "A", ]
  expect_true(all(ab$n == 64L))
})

test_that("planted pair is detected under per-program missingness (20% per program)", {
  fit <- .na_by_program(.cmp_toy(), n_miss = 16L, seed = 12)
  res <- suppressWarnings(find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none", seed = 12))
  expect_true(all(res$cmps$node_p < 0.05))
  expect_setequal(paste(res$cmps$cell_type, res$cmps$program_id), c("A A_1", "B B_2", "C C_3"))
  # every edge n is below the nominal 80 and at least 20 (min_shared)
  expect_true(all(res$edges$n < 80L))
  expect_true(all(res$edges$n >= 20L))
})
