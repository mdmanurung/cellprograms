# CP-2: assess_cmp_stability — synchronized donor resampling for CMPs.

## Toy fit: 3 cell types x 3 programs over 60 donors; program 1 of every type
## shares a latent donor driver (a planted CMP).
.make_cmpstab_fit <- function(N = 60, C = 3, K = 3, seed = 7, plant = TRUE) {
  set.seed(seed)
  ids <- paste0("d", seq_len(N))
  cts <- paste0("CT", seq_len(C))
  z <- if (plant) rnorm(N) else NULL
  scores <- lapply(cts, function(ct) {
    S <- matrix(rnorm(N * K), N, K, dimnames = list(ids, paste0(ct, "_P", seq_len(K))))
    if (plant) S[, 1] <- z + rnorm(N, sd = 0.3)
    S
  })
  names(scores) <- cts
  structure(list(
    scores = scores,
    loadings = lapply(cts, function(ct) matrix(0, 5, K)),
    cell_types = cts, observation_ids = ids,
    sample_metadata = data.frame(observation_id = ids, stringsAsFactors = FALSE)),
    class = c("cell_program_fit", "cell_program_data", "list"))
}

.cmpstab_args <- list(n_perm = 99, ambient = "none", min_shared = 10, min_donors = 10)

test_that("planted CMP is recovered under synchronized subsampling", {
  fit <- .make_cmpstab_fit()
  st <- assess_cmp_stability(fit, n_boot = 10, subsample_frac = 0.8,
                             cmp_args = .cmpstab_args, seed = 42, progress = FALSE)
  expect_s3_class(st, "cmp_stability")
  expect_equal(st$status, "ok")
  expect_equal(nrow(st$replicates), 10)
  expect_true(all(st$replicates$status == "ok"))
  # the planted 3-type CMP is in the reference and recovers frequently
  planted <- st$summary[vapply(strsplit(st$summary$members, ";"), length, 1) == 3, ]
  expect_gte(nrow(planted), 1)
  expect_gt(planted$recovery_freq[1], 0.6)
  expect_equal(planted$valid_refits[1], 10)
  expect_equal(planted$failed_refits[1], 0)
  expect_true(is.finite(planted$mean_jaccard[1]))
  # member-level table carries recovery and sign concordance
  expect_equal(nrow(st$members[st$members$cmp == planted$cmp[1], ]), 3)
  expect_true(all(st$cell_type_survival$survival_freq == 1))
})

test_that("no reference CMPs -> ineligibility reported, not zero stability", {
  fit <- .make_cmpstab_fit(plant = FALSE)
  expect_message(
    st <- assess_cmp_stability(fit, n_boot = 5, cmp_args = .cmpstab_args,
                               seed = 1, progress = FALSE),
    "no CMPs")
  expect_equal(st$status, "no_reference_cmps")
  expect_equal(nrow(st$summary), 0)
  expect_equal(nrow(st$replicates), 0)  # no replicates spent
})

test_that("donor_col with visits: cluster integrity, no pseudo-replication", {
  set.seed(11)
  donors <- paste0("d", 1:30)
  obs <- c(paste0(donors, "_v1"), paste0(donors, "_v2"))
  md <- data.frame(observation_id = obs,
                   donor = rep(donors, 2),
                   visit = rep(c("v1", "v2"), each = 30),
                   stringsAsFactors = FALSE)
  N <- length(obs)
  z <- rnorm(30)
  zobs <- z[match(md$donor, donors)]
  scores <- lapply(paste0("CT", 1:3), function(ct) {
    S <- matrix(rnorm(N * 2), N, 2, dimnames = list(obs, paste0(ct, "_P", 1:2)))
    S[, 1] <- zobs + rnorm(N, sd = 0.3)
    S
  })
  names(scores) <- paste0("CT", 1:3)
  fit <- structure(list(scores = scores,
                        loadings = lapply(scores, function(S) matrix(0, 5, 2)),
                        cell_types = names(scores), observation_ids = obs,
                        sample_metadata = md),
                   class = c("cell_program_fit", "cell_program_data", "list"))
  # visits present but visit_col missing -> explicit error naming the estimand
  expect_error(assess_cmp_stability(fit, donor_col = "donor", n_boot = 2,
                                    cmp_args = .cmpstab_args, progress = FALSE),
               "visit_col")
  # donor map + draw helpers: drawn donors keep ALL visits, ids stay unique
  map <- cellprograms:::.cmpstab_donor_map(fit, "donor", "visit", NULL)
  expect_equal(nrow(map), 60)
  expect_equal(length(unique(map$donor)), 30)
  set.seed(5)
  drawn <- cellprograms:::.cmpstab_draw(map$donor, 0.5, map$strata)
  expect_equal(length(drawn), 15)
  expect_equal(length(unique(drawn)), 15)  # without replacement
  obs_keep <- map$observation_id[map$donor %in% drawn]
  expect_equal(length(obs_keep), 30)       # 15 donors x 2 visits
  expect_equal(length(unique(obs_keep)), 30)  # no pseudo-replication
  fit_b <- cellprograms:::.cmpstab_build_fit(fit, obs_keep)
  expect_equal(length(fit_b$observation_ids), 30)
  expect_true(all(vapply(fit_b$scores, function(S) length(unique(rownames(S))) == nrow(S), TRUE)))
  # end-to-end with visits
  st <- assess_cmp_stability(fit, donor_col = "donor", visit_col = "visit",
                             n_boot = 5, subsample_frac = 0.8,
                             cmp_args = .cmpstab_args, seed = 3, progress = FALSE)
  expect_equal(st$status, "ok")
  expect_true(all(st$replicates$n_donors_drawn == 24))  # 0.8 * 30 donors
  expect_true(all(st$replicates$n_obs == 48))           # 24 donors x 2 visits
})

test_that("donor-level strata: draw is proportional within strata", {
  fit <- .make_cmpstab_fit()
  fit$sample_metadata$site <- rep(c("s1", "s2"), each = 30)
  map <- cellprograms:::.cmpstab_donor_map(fit, NULL, NULL, NULL)
  md <- fit$sample_metadata
  map$strata <- md$site[match(map$observation_id, md$observation_id)]
  set.seed(9)
  drawn <- cellprograms:::.cmpstab_draw(map$donor, 0.5, map$strata)
  tab <- table(md$site[match(drawn, md$observation_id)])
  expect_equal(as.integer(tab[["s1"]]), 15)
  expect_equal(as.integer(tab[["s2"]]), 15)
})

test_that("donor_strata must be constant within donor", {
  fit <- .make_cmpstab_fit()
  fit$sample_metadata$donor <- rep(paste0("d", 1:30), 2)[1:60]
  fit$sample_metadata$visit <- rep(c("v1", "v2"), each = 30)
  fit$sample_metadata$batch <- c(rep("x", 30), rep("y", 30))  # varies within each donor
  expect_error(cellprograms:::.cmpstab_donor_map(fit, "donor", "visit", "batch"),
               "constant within donor")
})

test_that("replicate failures are ledgered; rates stay NA on all-fail", {
  fit <- .make_cmpstab_fit()
  fake_ref <- structure(
    list(cmps = data.frame(cmp = 1L, cell_type = c("CT1", "CT2"),
                           program_id = c("CT1_P1", "CT2_P1"),
                           node_p = c(0.01, 0.01), sign = c(1, 1),
                           ambient_flag = c(NA, NA), stringsAsFactors = FALSE)),
    class = c("cmp_result", "list"))
  local_mocked_bindings(
    find_coordinated_programs = function(fit, ...) {
      if (nrow(fit$scores[[1]]) < 60) stop("mock replicate failure")
      fake_ref
    })
  st <- assess_cmp_stability(fit, n_boot = 4, subsample_frac = 0.8,
                             cmp_args = .cmpstab_args, seed = 2, progress = FALSE)
  expect_equal(st$status, "ok")
  expect_equal(nrow(st$replicates), 4)
  expect_true(all(st$replicates$status == "error"))
  expect_true(all(grepl("mock replicate failure", st$replicates$error_message)))
  expect_equal(st$summary$valid_refits[1], 0)
  expect_equal(st$summary$failed_refits[1], 4)
  expect_true(is.na(st$summary$recovery_freq[1]))  # ineligibility, not zero
})

test_that("incomplete donor x cell-type coverage: identical draw, honest eligibility", {
  fit <- .make_cmpstab_fit()
  # CT3 is missing rows for 15 donors entirely (no rows to synthesize)
  missing_ct3 <- paste0("d", 1:15)
  fit$scores$CT3 <- fit$scores$CT3[!rownames(fit$scores$CT3) %in% missing_ct3, , drop = FALSE]
  st <- assess_cmp_stability(fit, n_boot = 6, subsample_frac = 0.5,
                             cmp_args = .cmpstab_args, seed = 13, progress = FALSE)
  expect_equal(st$status, "ok")
  ok <- st$replicates[st$replicates$status == "ok", ]
  expect_equal(nrow(ok), 6)
  # every replicate drew 30 donors (0.5 * 60) exactly once
  expect_true(all(ok$n_donors_drawn == 30))
  # per-type eligible rows are ledgered and reflect incomplete coverage:
  # CT1/CT2 carry all 30 drawn donors; CT3 carries only the drawn donors
  # outside the 15 it lacks (never synthesized)
  elig <- strsplit(ok$eligible_rows, ";")
  ct1_n <- as.integer(sub("CT1:", "", sapply(elig, `[`, 1)))
  ct3_n <- as.integer(sub("CT3:", "", sapply(elig, `[`, 3)))
  expect_true(all(ct1_n == 30))
  expect_true(all(ct3_n <= 30) && all(ct3_n >= 30 - 15))
  # donor set applied to every type is identical by construction: replicate
  # fits subset one shared obs_keep vector (verified via helper directly)
  map <- cellprograms:::.cmpstab_donor_map(fit, NULL, NULL, NULL)
  set.seed(99)
  drawn <- cellprograms:::.cmpstab_draw(map$donor, 0.5, map$strata)
  obs_keep <- map$observation_id[map$donor %in% drawn]
  fit_b <- cellprograms:::.cmpstab_build_fit(fit, obs_keep)
  expect_true(all(rownames(fit_b$scores$CT1) %in% obs_keep))
  expect_true(all(rownames(fit_b$scores$CT3) %in% obs_keep))
  expect_equal(sort(rownames(fit_b$scores$CT1)), sort(obs_keep))
  expect_equal(sort(rownames(fit_b$scores$CT3)), sort(setdiff(obs_keep, missing_ct3)))
})

test_that("seeded runs are reproducible", {
  fit <- .make_cmpstab_fit()
  a <- assess_cmp_stability(fit, n_boot = 5, cmp_args = .cmpstab_args, seed = 77, progress = FALSE)
  b <- assess_cmp_stability(fit, n_boot = 5, cmp_args = .cmpstab_args, seed = 77, progress = FALSE)
  expect_identical(a$summary$recovery_freq, b$summary$recovery_freq)
  expect_identical(a$replicates$n_donors_drawn, b$replicates$n_donors_drawn)
})

test_that("bad inputs fail fast", {
  fit <- .make_cmpstab_fit()
  expect_error(assess_cmp_stability(fit, donor_col = "nope", cmp_args = .cmpstab_args),
               "donor_col 'nope'")
  expect_error(assess_cmp_stability(fit, subsample_frac = 1.5, cmp_args = .cmpstab_args),
               "subsample_frac")
  expect_error(assess_cmp_stability(fit, cmp_args = list(fit = fit)), "cmp_args")
})
