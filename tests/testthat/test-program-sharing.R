# Toy canonicalized fit with planted loading structure and posterior SDs;
# no flashier needed (loadings/SDs planted directly, .cmp_toy style).
# A_1 = B_1 (identical shared program), A_2 private to A, C_1 private to C.
# null_rho plants correlated null-gene loading noise (shared estimation noise).
.sharing_toy <- function(g = 300, seed = 1, null_rho = 0) {
  set.seed(seed)
  genes <- sprintf("G%03d", seq_len(g))
  cts <- c("A", "B", "C")
  # Effects bounded away from zero (|e| in [1.5, 3], z >= 5): the tests assert
  # sharing CLASSIFICATION, so planted signal must not have an undetectable
  # weak tail (with rnorm(sd = 2) ~25% of genes have |z| < 2 and mash --
  # correctly, given ~87% of genes are truly null -- leaves them uncalled).
  strong <- function() runif(40, 1.5, 3) * sample(c(-1, 1), 40, replace = TRUE)
  eff <- list(
    A = list(c(strong(), rep(0, g - 40)),                            # A_1: shared
             c(rep(0, 60), strong(), rep(0, g - 100))),              # A_2: private
    B = list(NULL),
    C = list(c(rep(0, 120), strong(), rep(0, g - 160)))              # C_1: private
  )
  eff$B[[1]] <- eff$A[[1]]                                         # identical shared program
  loadings <- fits <- canon <- list()
  common <- rnorm(g)   # shared across cell types: plants CROSS-CELL-TYPE null correlation
  for (ct in cts) {
    K <- length(eff[[ct]])
    W <- matrix(rnorm(g * K, sd = 0.05), g, K,
                dimnames = list(genes, paste0(ct, "_", seq_len(K))))
    for (k in seq_len(K)) {
      W[, k] <- W[, k] + eff[[ct]][[k]] + sqrt(null_rho) * 0.05 * common
    }
    loadings[[ct]] <- W
    fits[[ct]] <- list(
      flash = list(F_pm = W, F_psd = matrix(0.3, g, K, dimnames = dimnames(W))),
      genes = genes
    )
    canon[[ct]] <- list(scale_multiplier = rep(1, K), sign_multiplier = rep(1, K))
  }
  structure(list(loadings = loadings, fits = fits, canonicalization = canon,
                 cell_types = cts),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

test_that("matching recovers a planted shared program and rejects unrelated ones", {
  fit <- .sharing_toy()
  m <- .match_programs_loadings(fit, n_perm = 199, seed = 1)
  expect_true(all(c("ct_a", "prog_a", "ct_b", "prog_b", "cor", "p", "q") %in% names(m)))
  hit <- m[m$ct_a == "A" & m$prog_a == "A_1" & m$ct_b == "B" & m$prog_b == "B_1", ]
  expect_equal(nrow(hit), 1L)
  expect_gt(abs(hit$cor), 0.9)
  expect_lt(hit$q, 0.05)
  # no other program pair (all unrelated) is called a match
  expect_false(any(m$q[-match(rownames(hit), rownames(m))] < 0.05))
  # grouping puts A_1 + B_1 together and nothing else
  gr <- .sharing_groups(m, alpha = 0.05)
  expect_equal(length(gr), 1L)
  expect_setequal(gr[[1]]$program_id, c("A_1", "B_1"))
})

test_that("V estimation returns a valid shrunk correlation matrix", {
  fit <- .sharing_toy(null_rho = 0.5)   # planted null-noise correlation ~ rho/(1+rho) = 1/3
  grp <- data.frame(cell_type = c("A", "B"), program_id = c("A_1", "B_1"))
  V <- .estimate_V(fit, grp)
  expect_true(isSymmetric(V))
  expect_equal(unname(diag(V)), c(1, 1))
  expect_equal(dimnames(V), list(c("A", "B"), c("A", "B")))
  expect_gt(attr(V, "n_null"), 200L)
  expect_gt(V[1, 2], 0.1)               # recovers the planted null correlation
  expect_lt(V[1, 2], 0.6)
  expect_true(all(eigen(V)$values > 0))
  # too few null genes -> identity + warning
  fit2 <- fit
  fit2$loadings$A[, 1] <- fit2$loadings$A[, 1] * 100   # no null genes left
  fit2$loadings$B[, 1] <- fit2$loadings$B[, 1] * 100
  expect_warning(V2 <- .estimate_V(fit2, grp), "identity")
  expect_equal(V2, diag(2), ignore_attr = TRUE)
})

test_that("program_sharing errors gracefully without mashr or without posterior SDs", {
  fit <- .sharing_toy()
  testthat::local_mocked_bindings(.has_mashr = function() FALSE,
                                  .package = "cellprograms")
  expect_error(program_sharing(fit), "mashr")
})

test_that("program_sharing tells users to refit when F_psd is missing", {
  skip_if_not_installed("mashr")
  fit <- .sharing_toy()
  fit$fits$A$flash$F_psd <- NULL
  expect_error(program_sharing(fit), "Refit")
})

test_that("degenerate zero-SD entries: floored when loading is zero, dropped when inconsistent", {
  skip_if_not_installed("mashr")
  # zero loading + zero SD (program does not use the gene): floored, gene kept
  fit <- .sharing_toy()
  g5 <- sprintf("G%03d", 41:45)              # null genes for the shared A_1/B_1 pair
  fit$loadings$A[g5, "A_1"] <- 0; fit$loadings$B[g5, "B_1"] <- 0
  fit$fits$A$flash$F_psd[g5, ] <- 0; fit$fits$B$flash$F_psd[g5, ] <- 0
  expect_warning(res <- program_sharing(fit, n_perm = 199, seed = 1), "floor")
  expect_equal(nrow(res$results$group1$posterior_mean), 300L)
  expect_equal(res$results$group1$classification, "shared")
  # nonzero loading + zero SD (inconsistent posterior): dropped
  fit2 <- .sharing_toy()
  fit2$fits$A$flash$F_psd[1:5, ] <- 0
  fit2$fits$B$flash$F_psd[1:5, ] <- 0
  expect_warning(res2 <- program_sharing(fit2, n_perm = 199, seed = 1), "Dropping")
  expect_equal(nrow(res2$results$group1$posterior_mean), 295L)
  # z_thresh pass-through: wider null window finds null genes where z_thresh = 1 finds none.
  # Scale SDs so null-gene z ~ N(0, 1); at z_thresh = 1 only ~48% of null genes qualify
  # (both conditions), at 3 nearly all.
  fit3 <- .sharing_toy(null_rho = 0.5)
  fit3$fits$A$flash$F_psd <- fit3$fits$A$flash$F_psd * 0.2
  fit3$fits$B$flash$F_psd <- fit3$fits$B$flash$F_psd * 0.2
  grp2 <- data.frame(cell_type = c("A", "B"), program_id = c("A_1", "B_1"))
  V1 <- .estimate_V(fit3, grp2, z_thresh = 1)
  V2 <- .estimate_V(fit3, grp2, z_thresh = 3)
  expect_gt(attr(V2, "n_null"), attr(V1, "n_null"))
})

test_that("numerically degenerate programs are excluded from matching and sharing", {
  skip_if_not_installed("mashr")
  # crush A_2 into numerical noise (15 orders below A_1): matching skips it
  fit <- .sharing_toy()
  fit$loadings$A[, "A_2"] <- fit$loadings$A[, "A_2"] * 1e-16
  expect_warning(mt <- .match_programs_loadings(fit), "degenerate")
  expect_false(any(mt$prog_a == "A_2" | mt$prog_b == "A_2"))
  expect_true(any(mt$prog_a == "A_1" | mt$prog_b == "A_1"))
  # an explicitly supplied degenerate member is dropped from the matrices
  grp <- data.frame(cell_type = c("A", "B", "C"), program_id = c("A_2", "B_1", "C_1"))
  expect_warning(mats <- .sharing_matrices(fit, grp), "degenerate")
  expect_equal(ncol(mats$B), 2L)
  expect_equal(colnames(mats$B), c("B", "C"))
  # a group reduced to < 2 non-degenerate members is an error
  grp2 <- data.frame(cell_type = c("A", "B"), program_id = c("A_2", "B_1"))
  expect_error(.sharing_matrices(fit, grp2), "fewer than 2")
  # a real-scale program with a collapsed posterior-SD column (fixed factor)
  # is degenerate too
  fit4 <- .sharing_toy()
  fit4$fits$A$flash$F_psd[, "A_2"] <- 1e-20
  expect_warning(mt4 <- .match_programs_loadings(fit4), "degenerate")
  expect_false(any(mt4$prog_a == "A_2" | mt4$prog_b == "A_2"))
})

test_that("end-to-end: planted shared program is shared, private programs are private", {
  skip_if_not_installed("mashr")
  fit <- .sharing_toy()
  res <- program_sharing(fit, n_perm = 199, seed = 1)
  expect_s3_class(res, "program_sharing")
  s <- res$summary
  expect_equal(s$classification[s$program_id == "A_1"], "shared")
  expect_equal(s$classification[s$program_id == "B_1"], "shared")
  expect_equal(s$classification[s$program_id == "A_2"], "private")
  expect_equal(s$classification[s$program_id == "C_1"], "private")
  expect_equal(s$evidence[s$program_id == "C_1"], "no_significant_match")
  r <- res$results$group1
  expect_equal(dim(r$posterior_mean), c(300L, 2L))
  expect_equal(dim(r$lfsr), c(300L, 2L))
  expect_gt(r$p_all_shared, 0.5)
  expect_gt(sum(r$pattern_weights), 0.99)
  expect_true("equal_effects" %in% names(r$pattern_weights))
  # forced match of two unrelated private programs -> mash calls it private
  forced <- list(data.frame(cell_type = c("A", "C"), program_id = c("A_2", "C_1")))
  res2 <- program_sharing(fit, matches = forced)
  expect_equal(res2$results$group1$classification, "private")
  expect_gte(res2$results$group1$frac_top_private, 0.8)   # ~5% of null-condition lfsr calls are false by construction
})

test_that("fit_celltype_programs stores loading posterior SDs (both modes)", {
  skip_if_not_installed("flashier")
  skip_if_not_installed("ebnm")
  fit <- fit_celltype_programs(.toy(), max_factors = 2, seed = 1)
  fl <- fit$fits$A$flash
  expect_false(is.null(fl$F_psd))
  expect_equal(dim(fl$F_psd), dim(fl$F_pm))
  expect_equal(dim(fl$L_psd), dim(fl$L_pm))
  fit_f <- fit_celltype_programs(.toy(), max_factors = 2, covariates = "batch",
                                 covariate_mode = "fixed", seed = 1)
  flf <- fit_f$fits$A$flash
  expect_false(is.null(flf$F_psd))
  expect_equal(ncol(flf$F_psd), ncol(flf$F_pm))
})

test_that("program_sharing with zero sharing groups returns an all-private summary", {
  skip_if_not_installed("mashr")
  # Cell types A and C share no planted program (A_1/A_2 vs C_1 unrelated), so
  # matching finds zero significant pairs. Regression: names(results) <-
  # paste0("group", seq_along(results)) crashed on the empty results list
  # (paste0 treats integer(0) as "" -> length-1 names on a length-0 list).
  fit <- .sharing_toy()
  fit$loadings <- fit$loadings[c("A", "C")]
  fit$fits <- fit$fits[c("A", "C")]
  fit$canonicalization <- fit$canonicalization[c("A", "C")]
  fit$cell_types <- c("A", "C")
  res <- program_sharing(fit, n_perm = 199, seed = 1)
  expect_length(res$groups, 0L)
  expect_length(res$results, 0L)
  expect_true(all(res$summary$classification == "private"))
  expect_true(all(res$summary$evidence == "no_significant_match"))
  expect_equal(nrow(res$summary), 3L)   # A_1, A_2, C_1
})
