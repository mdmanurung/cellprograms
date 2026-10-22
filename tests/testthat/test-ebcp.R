## EB-CP tensor factorization: recovery, trap resistance, fit-wrapper path.
##
## NB: simulate_pseudobulk() plants DISJOINT gene sets per cell type for its
## shared-activity programs — the regime where eb_cp's shared gene mode is
## structurally diluted (documented in benchmarks/coordination/PREREG.md §7).
## The positive controls here therefore plant shared activity AND shared genes
## (as helper-toy::.toy() does), matching the regime the C1 gate test validated.

## hand-planted act+genes matrices: A and B share activity z AND loading w;
## C has an independent activity on the same genes (private)
make_mats <- function(seed, signal = 1.2) {
  set.seed(seed)
  G <- 200; donors <- paste0("d", seq_len(40))
  z <- stats::rnorm(40)
  w <- c(rep(0, 150), stats::rnorm(50))
  mats <- list()
  for (ct in c("A", "B", "C")) {
    Y <- matrix(stats::rnorm(40 * G), 40, G,
                dimnames = list(donors, paste0("g", seq_len(G))))
    if (ct %in% c("A", "B")) Y <- Y + signal * z %*% t(w)
    if (ct == "C") Y <- Y + signal * stats::rnorm(40) %*% t(w)
    mats[[ct]] <- Y
  }
  mats
}

test_that("fit_eb_cp recovers a shared activity+genes program", {
  cp <- fit_eb_cp(make_mats(11), Kmax = 8L, seed = 1L)
  expect_s3_class(cp, "eb_cp_fit")
  expect_gt(cp$n_components, 0L)
  pr <- eb_cp_programs(cp)
  coord <- pr$programs[pr$programs$is_coordination, ]
  expect_true(any(grepl("\\bA\\b", coord$cell_type) & grepl("\\bB\\b", coord$cell_type)),
              label = "a coordination program with both planted members A and B")
})

test_that("fit_eb_cp does not call a coordination program on correlated privates", {
  ## A and B get correlated (r = 0.7) donor activities on DISJOINT genes: the
  ## S5 trap. Truth: distinct programs; a call containing both is a failure.
  set.seed(11)
  G <- 200; donors <- paste0("d", seq_len(40))
  zA <- stats::rnorm(40)
  zB <- 0.7 * zA + sqrt(1 - 0.49) * stats::rnorm(40)
  wA <- c(rep(0, 100), stats::rnorm(100))
  wB <- c(stats::rnorm(100), rep(0, 100))
  mats <- list(
    A = matrix(stats::rnorm(40 * G), 40, G, dimnames = list(donors, paste0("g", seq_len(G)))) +
      1.2 * zA %*% t(wA),
    B = matrix(stats::rnorm(40 * G), 40, G, dimnames = list(donors, paste0("g", seq_len(G)))) +
      1.2 * zB %*% t(wB)
  )
  cp <- fit_eb_cp(mats, Kmax = 8L, seed = 1L)
  pr <- eb_cp_programs(cp)
  trap <- pr$programs$is_coordination &
    grepl("\\bA\\b", pr$programs$cell_type) & grepl("\\bB\\b", pr$programs$cell_type)
  expect_false(any(trap), label = "no coordination call over the correlated-private pair")
})

test_that("eb_cp_from_fit derives strata from the fit's covariates", {
  x <- .toy()   # A + B share a latent program; batch/age nuisance; metadata attached
  f <- canonicalize_programs(fit_celltype_programs(
    x, max_factors = 6L, covariates = "batch", covariate_mode = "residualize",
    seed = 3L))
  cp <- eb_cp_from_fit(f, Kmax = 8L, seed = 1L)
  expect_s3_class(cp, "eb_cp_fit")
  n_ids <- length(unique(unlist(lapply(x$pseudobulk, rownames))))
  expect_length(cp$strata_used, n_ids)             # derived, one per donor
  expect_identical(sort(unique(cp$strata_used)), c("a", "b"))
  pr <- eb_cp_programs(cp)
  expect_true(any(pr$programs$is_coordination),
              label = "the A+B shared program is a coordination call")
})

test_that("fit_eb_cp validates inputs", {
  m <- matrix(0, 5, 5)
  expect_error(fit_eb_cp(list(m, m), seed = 1), "named list")
  m3 <- matrix(stats::rnorm(60), 20, 3, dimnames = list(paste0("d", 1:20), paste0("g", 1:3)))
  expect_error(fit_eb_cp(list(A = m3, B = m3), seed = 1), "shared")
})
