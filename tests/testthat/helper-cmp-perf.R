# Fixtures for the find_coordinated_programs() regression guard
# (test-cmp-perf-regression.R). Scores are planted directly; no flashier needed.
# The golden outputs in fixtures/cmp_perf_golden.rds were captured from the
# permutation loop as it was before the fast path (release/0.3.0 @ 0914144).

.perf_fit <- function(n = 40, k = 3, g = 100, seed = 2, cts = c("A", "B", "C", "D")) {
  set.seed(seed)
  n_pat <- n / 2
  patient <- rep(sprintf("p%02d", seq_len(n_pat)), each = 2)
  visit <- rep(c("T1", "T2"), n_pat)
  ids <- paste(patient, visit, sep = "_")
  site <- rep(c("s1", "s2"), length.out = n_pat)[match(patient, unique(patient))]
  ctx <- rep(c(0, 1), length.out = n)
  zp <- rnorm(n_pat)
  z <- zp[match(patient, unique(patient))] + rnorm(n, sd = 0.3)
  scores <- loadings <- pseudobulk <- list()
  for (ct in cts) {
    S <- matrix(rnorm(n * k), n, k, dimnames = list(ids, paste0(ct, "_", seq_len(k))))
    j <- match(ct, cts)
    # A_1, B_2 and C_3 share one driver (C with the opposite sign); D is uncoupled
    if (j <= 3) S[, j] <- z * (if (j == 3) -1 else 1) + rnorm(n, sd = 0.6)
    W <- matrix(rnorm(g * k, sd = 0.1), g, k,
                dimnames = list(sprintf("G%03d", seq_len(g)), colnames(S)))
    W[(1:15) + 15 * (j - 1), ] <- W[(1:15) + 15 * (j - 1), ] + 1
    scores[[ct]] <- S
    loadings[[ct]] <- W
    pseudobulk[[ct]] <- matrix(rnorm(n * g, 5), n, g, dimnames = list(ids, rownames(W)))
  }
  structure(list(scores = scores, loadings = loadings, cell_types = cts, observation_ids = ids,
                 pseudobulk = pseudobulk,
                 sample_metadata = data.frame(observation_id = ids, patient = patient, visit = visit,
                                              site = site, ctx = ctx, stringsAsFactors = FALSE)),
            class = c("cell_program_fit", "cell_program_data", "list"))
}

.perf_cases <- function() {
  base <- .perf_fit()
  na_fit <- base
  na_fit$scores$A[1:6, 1] <- NA
  na_fit$scores$B[5:12, 2] <- NA
  sparse_fit <- base
  sparse_fit$scores$C[26:40, ] <- NA      # C observed in 25 donors only
  const_fit <- base
  const_fit$scores$A[, 3] <- 0            # constant program: constant_residual
  list(
    complete = list(fit = base, args = list(strata = "site", ambient = "none")),
    ambient_flag = list(fit = base, args = list(strata = "site", ambient = "flag")),
    within_context = list(fit = base, args = list(strata = "site", context = "ctx", mode = "within",
                                                  ambient = "none")),
    clusters = list(fit = base, args = list(strata = "site", clusters = "patient", visits = "visit",
                                            ambient = "none")),
    na_donors = list(fit = na_fit, args = list(strata = "site", ambient = "none")),
    sparse_celltype = list(fit = sparse_fit, args = list(strata = "site", ambient = "none",
                                                         min_shared = 30L, min_donors = 20L)),
    constant_program = list(fit = const_fit, args = list(strata = "site", ambient = "none"))
  )
}

.perf_run <- function(case, n_perm = 99L, seed = 1L) {
  suppressWarnings(do.call(find_coordinated_programs,
                           c(list(fit = case$fit, n_perm = n_perm, seed = seed), case$args)))
}

.perf_pin <- function(res) {
  list(nodes = res$nodes, cmps = res$cmps, edges = res$edges, edge_support = res$edge_support)
}

.perf_stability_pin <- function() {
  st <- suppressWarnings(assess_cmp_stability(.perf_fit(), n_boot = 3, subsample_frac = 0.8,
                                              cmp_args = list(strata = "site", ambient = "none",
                                                              n_perm = 99L),
                                              seed = 1L, progress = FALSE))
  # elapsed_seconds is wall time, not a result
  list(summary = st$summary, members = st$members,
       replicates = st$replicates[setdiff(names(st$replicates), "elapsed_seconds")])
}
