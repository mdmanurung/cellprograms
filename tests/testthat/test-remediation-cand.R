# Regression tests for the 2026-10-18 remediation queue:
# CAND-01, CAND-05, CAND-07 (docs), D-2a, CAND-03+D-3b, CAND-12.
# Fix scopes registered in reports/audit/TASKS.md before any R/ change.

# ---- CAND-03: non-finite contract errors --------------------------------
test_that("as_cell_program_data rejects Inf/NaN, not just NA (CAND-03)", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 60, program_size = 10, seed = 1)
  pb <- sim$pseudobulk
  ct1 <- names(pb)[1]
  pb_inf <- pb; pb_inf[[ct1]][1, 1] <- Inf
  expect_error(as_cell_program_data(pb_inf), "non-finite")
  pb_nan <- pb; pb_nan[[ct1]][1, 1] <- NaN
  expect_error(as_cell_program_data(pb_nan), "non-finite")
  pb_na <- pb; pb_na[[ct1]][1, 1] <- NA
  expect_error(as_cell_program_data(pb_na), "non-finite")
})

# ---- D-3b: .greedy_match survives non-finite similarities ---------------
test_that(".greedy_match demotes non-finite similarities instead of crashing (D-3b)", {
  S <- matrix(c(1, NaN, 0.5, 0.9, Inf, 0.2, 0.3, 0.4, NA), 3, 3)
  expect_warning(m <- .greedy_match(S), "non-finite")
  expect_true(all(is.finite(m$sim)))
  expect_equal(m$sim[1], 1)  # global finite maximum matched first
  expect_identical(nrow(m), 3L)
  S2 <- matrix(NaN, 2, 2)
  expect_warning(m2 <- .greedy_match(S2), "non-finite")
  expect_identical(nrow(m2), 0L)  # nothing matchable, no crash
})

# ---- CAND-01: stability failure accounting ------------------------------
test_that("stability refit failures are counted and excluded from the denominator (CAND-01)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  orig <- fit_celltype_programs
  calls <- 0L
  ## fail every second refit deterministically: 2 of 4 fail
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    calls <<- calls + 1L
    if (calls %% 2L == 0L) stop("injected refit failure")
    orig(...)
  })
  expect_warning(fit_s <- assess_program_stability(fit, n_boot = 4, seed = 1, progress = FALSE),
                 "2 of 4 subsample refits failed")
  st <- fit_s$stability$per_program
  expect_identical(fit_s$stability$n_failed, 2L)
  expect_identical(fit_s$stability$n_success, 2L)
  ## the registered contract: recovery_freq == n_hit / n_success exactly
  expect_equal(st$recovery_freq, st$n_hit / fit_s$stability$n_success)
  ## a program recovered in both successful refits reads 1.0, not 0.5
  expect_equal(max(st$recovery_freq), 1)
})

test_that("all-fail stability is NA with a distinct warning, not 0 (CAND-01)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) stop("injected"))
  expect_warning(fit_s <- assess_program_stability(fit, n_boot = 2, seed = 1, progress = FALSE),
                 "all 2 subsample refits failed")
  expect_true(all(is.na(fit_s$stability$per_program$recovery_freq)))
  expect_identical(fit_s$stability$n_success, 0L)
})

# ---- CAND-05: chain orientation propagation -----------------------------
.cmp_chain_toy <- function(n = 80, seed = 7) {
  ## A_1 ~ +z1 (strong), B_1 ~ +z1 + z2 (strong link to both), C_1 ~ -z2
  ## (strong ANTI link to B_1); the direct A_1-C_1 edge is noise (|T| < 3)
  ## with a random sign. The CMP forms via the A-B and B-C linking edges.
  set.seed(seed)
  ids <- sprintf("d%03d", seq_len(n))
  site <- rep(c("s1", "s2"), length.out = n)
  z1 <- rnorm(n); z2 <- rnorm(n)
  genes <- sprintf("G%03d", seq_len(200))
  mk_S <- function(ct, col) {
    S <- matrix(rnorm(n * 3, sd = 0.3), n, 3, dimnames = list(ids, paste0(ct, "_", 1:3)))
    S[, col] <- switch(col, z1 + rnorm(n, sd = 0.3),
                       z1 + z2 + rnorm(n, sd = 0.3),
                       -z2 + rnorm(n, sd = 0.3))
    S
  }
  scores <- loadings <- pb <- list()
  for (ct in c("A", "B", "C")) {
    S <- mk_S(ct, match(ct, c("A", "B", "C")))
    W <- matrix(rnorm(200 * 3, sd = 0.1), 200, 3, dimnames = list(genes, colnames(S)))
    scores[[ct]] <- S; loadings[[ct]] <- W
    pb[[ct]] <- matrix(rnorm(n * 200, 5), n, 200, dimnames = list(ids, genes))
  }
  structure(list(scores = scores, loadings = loadings, cell_types = c("A", "B", "C"),
                 observation_ids = ids, pseudobulk = pb,
                 sample_metadata = data.frame(observation_id = ids, site = site,
                                              stringsAsFactors = FALSE)),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

test_that("chain-joined CMP members orient along linking edges, not noise edges (CAND-05)", {
  fit <- .cmp_chain_toy(seed = 7)
  res <- find_coordinated_programs(fit, strata = "site", n_perm = 99, ambient = "none", seed = 7)
  cmps <- res$cmps
  expect_true(all(c("A", "B", "C") %in% cmps$cell_type))  # one 3-member CMP formed
  s <- setNames(cmps$sign, paste0(cmps$cell_type, "|", cmps$program_id))
  sA <- s[startsWith(names(s), "A|")][1]
  sB <- s[startsWith(names(s), "B|")][1]
  sC <- s[startsWith(names(s), "C|")][1]
  ## path products: A-B linked positive, B-C linked ANTI. Root-invariant.
  expect_equal(unname(sB), unname(sA))   # product along A-B path = +
  expect_equal(unname(sC), unname(-sA))  # product along A-B-C path = -
})

# ---- D-2a/D-2c: matches validation + forced-pair warning ----------------
test_that("program_sharing validates matches references and shape (D-2a/D-2c)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  bad_ref <- data.frame(cell_type = c("A", "B"), program_id = c("A_1", "NOPE"))
  expect_error(program_sharing(fit, matches = bad_ref), "unknown program")
  bad_shape <- data.frame(foo = 1:3, bar = 2:4)
  expect_error(program_sharing(fit, matches = bad_shape), "matched set")
})

test_that("forced matches emit the D-2a interpretation warning (D-2a)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  skip_if_not_installed("mashr")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  forced <- data.frame(cell_type = c("A", "B"), program_id = c("A_1", "B_1"))
  expect_warning(res <- program_sharing(fit, matches = forced, n_top = 20), "user-forced")
  expect_true("classification" %in% names(res$summary))
})

# ---- CAND-12: spec storage, replay, conflict list -----------------------
test_that("frozen spec is stored on the fit and replayed by stability refits (CAND-12)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  spec <- freeze_preprocessing(.toy(), covariates = ~ batch, features = "variable",
                               n_variable_genes = 150)
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), preprocessing_spec = spec,
                                                     max_factors = 3, seed = 1))
  expect_equal(fit$preprocessing$frozen_spec, spec$digest)
  expect_false(is.null(fit$preprocessing$spec))
  seen <- list(); orig <- fit_celltype_programs
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...); seen[[length(seen) + 1L]] <<- names(a); orig(...)
  })
  assess_program_stability(fit, n_boot = 1, seed = 1, progress = FALSE)
  expect_length(seen, 1L)
  expect_true("preprocessing_spec" %in% seen[[1]])  # replay path
  expect_false("features" %in% seen[[1]])           # spec governs; no explicit conflicts
  expect_false("center" %in% seen[[1]])
})

test_that("feature_scale is in the MR-0 conflict list and passes through legacy refits (CAND-12)", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  spec <- freeze_preprocessing(.toy(), covariates = ~ batch, features = "variable",
                               n_variable_genes = 150)
  expect_error(fit_celltype_programs(.toy(), preprocessing_spec = spec,
                                     feature_scale = "log2", max_factors = 3, seed = 1),
               "frozen")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  seen <- list(); orig <- fit_celltype_programs
  testthat::local_mocked_bindings(fit_celltype_programs = function(...) {
    a <- list(...); seen[[length(seen) + 1L]] <<- a$feature_scale; orig(...)
  })
  assess_program_stability(fit, n_boot = 1, seed = 1, progress = FALSE)
  expect_identical(seen[[1]], fit$preprocessing$feature_scale)
})
