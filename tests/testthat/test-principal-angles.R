test_that("principal_angles on adjusted fits finds the shared program", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), covariates = c("batch", "age"),
                                                     max_factors = 3, seed = 1))
  pa <- principal_angles(fit, "A", "B", n_perm = 99, seed = 1)
  expect_gt(pa$cosines[1], 0.9)
  expect_lt(pa$p_values[1], 0.05)
  pl <- .expect_warning_messages(
    principal_angles(fit, "A", "B", space = "loadings", weight = "residual_sd", n_perm = 0),
    "Gene intersection for A vs B is 60 (< 200); angles may be unstable.")
  expect_s3_class(pl, "principal_angles")
})

test_that(".perm_within permutes only inside strata", {
  set.seed(1)
  st <- rep(c("x", "y", "z"), c(5, 6, 7))
  for (i in 1:20) {
    p <- .perm_within(length(st), st)
    expect_identical(st[p], st)
    expect_setequal(p, seq_along(st))
  }
})

test_that("strata is validated and accepted", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  ids <- intersect(rownames(fit$scores$A), rownames(fit$scores$B))
  st <- setNames(rep(c("u", "v"), length.out = length(ids)), ids)
  pa <- principal_angles(fit, "A", "B", n_perm = 19, seed = 1, strata = st)
  expect_s3_class(pa, "principal_angles")
  expect_error(principal_angles(fit, "A", "B", n_perm = 5, strata = unname(st)), "named vector")
  expect_error(principal_angles(fit, "A", "B", space = "loadings", n_perm = 5, strata = st), "scores space")
})

test_that("sharing_spectrum gives a named error when fewer than 2 cell types retain a program", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- canonicalize_programs(fit_celltype_programs(.toy(), max_factors = 3, seed = 1))
  fit1 <- fit
  fit1$loadings$B <- fit1$loadings$B[, 0L, drop = FALSE]
  expect_error(sharing_spectrum(fit1, n_perm = 19), "found 1")
  for (ct in c("A", "B")) fit$loadings[[ct]] <- fit$loadings[[ct]][, 0L, drop = FALSE]
  expect_error(sharing_spectrum(fit, n_perm = 19), "found 0")
})

# ---- repeated visits: cluster (patient) permutation -----------------------------------

test_that("default calls keep the RNG stream of the pre-change code at seed = 1", {
  skip_if_not_installed("flashier"); skip_if_not_installed("ebnm")
  fit <- .pa_golden_fit()
  md <- fit$sample_metadata
  cts <- fit$cell_types[1:2]
  # golden values captured from release/0.3.0 @ 166e40c, before this change
  a <- principal_angles(fit, cts[1], cts[2], n_perm = 99, seed = 1)
  expect_equal(c(a$p_pair, a$z_pair), c(0.01, 18.9515028894), tolerance = 1e-8)
  b <- principal_angles(fit, cts[1], cts[2], n_perm = 99, seed = 1,
                        strata = stats::setNames(md$site, md$observation_id))
  expect_equal(c(b$p_pair, b$z_pair), c(0.01, 18.9576037446), tolerance = 1e-8)
  l <- principal_angles(fit, cts[1], cts[2], space = "loadings", n_perm = 99, seed = 1)
  expect_equal(c(l$p_pair, l$z_pair), c(0.28, 0.0707608442685), tolerance = 1e-8)
  # seed = 1L is now the default
  expect_equal(principal_angles(fit, cts[1], cts[2], n_perm = 99)[c("p_pair", "z_pair")],
               a[c("p_pair", "z_pair")])
})

test_that("clusters and visits validate their design", {
  fit <- .cmp_visit_toy(n_pat = 12)
  pa <- function(...) principal_angles(fit, "A", "B", n_perm = 19, ...)
  expect_error(pa(clusters = "patient"), "supply `visits`")
  expect_error(pa(visits = "visit"), "requires `clusters`")
  expect_error(pa(clusters = "patient", visits = "visit", space = "loadings"), "scores space")
  expect_error(pa(strata = stats::setNames(fit$sample_metadata$visit_site, fit$observation_ids),
                  clusters = "patient", visits = "visit"), "constant within each cluster")
  expect_error(pa(clusters = "no_such_column", visits = "visit"), "not found")
})

test_that("one-observation clusters reproduce the row permutation exactly", {
  fit <- .cmp_visit_toy(n_pat = 20, visits = "T1")
  a <- principal_angles(fit, "A", "B", n_perm = 49, seed = 3)
  b <- principal_angles(fit, "A", "B", n_perm = 49, seed = 3, clusters = "observation_id")
  expect_equal(a$p_values, b$p_values)
  expect_equal(a$z_pair, b$z_pair)
})

test_that("cluster permutation moves whole patients, so observed angles are unchanged", {
  fit <- .cmp_visit_toy(n_pat = 20)
  a <- principal_angles(fit, "A", "B", n_perm = 49, seed = 3)
  b <- principal_angles(fit, "A", "B", n_perm = 49, seed = 3, clusters = "patient", visits = "visit")
  expect_equal(a$angles, b$angles)
  expect_false(isTRUE(all.equal(a$z_pair, b$z_pair)))   # a different null
})

test_that("repeated visits: cluster permutation is less anti-conservative than row permutation", {
  p <- vapply(1:40, function(s) {
    fit <- .cmp_visit_toy(n_pat = 30, seed = 700 + s, cts = c("A", "B"))
    c(row = principal_angles(fit, "A", "B", n_perm = 99, seed = s)$p_pair,
      cluster = principal_angles(fit, "A", "B", n_perm = 99, seed = s,
                                 clusters = "patient", visits = "visit")$p_pair)
  }, c(row = 0, cluster = 0))
  rej <- rowSums(p < 0.05)
  expect_lte(rej[["cluster"]], rej[["row"]])
  expect_gt(rej[["row"]], 0.05 * 40)          # row permutation rejects above its nominal rate
})

test_that("too few exchangeable clusters is a classed error and few distinct permutations warn", {
  few <- .cmp_visit_toy(n_pat = 1)
  err <- tryCatch(principal_angles(few, "A", "B", n_perm = 19, clusters = "patient", visits = "visit"),
                  error = identity)
  expect_s3_class(err, "cp_few_clusters")
  expect_match(conditionMessage(err), "fewer than 2")
  small <- .cmp_visit_toy(n_pat = 4)
  expect_warning(principal_angles(small, "A", "B", n_perm = 99, clusters = "patient", visits = "visit"),
                 "distinct cluster permutations")
})

test_that("sharing_spectrum forwards clusters and visits and skips pairs it cannot test", {
  fit <- .cmp_visit_toy(n_pat = 25)
  ss <- sharing_spectrum(fit, n_perm = 49, seed = 1, clusters = "patient", visits = "visit")
  expect_equal(nrow(ss$summary), 3L)
  expect_error(sharing_spectrum(fit, n_perm = 49, seed = 1, clusters = "patient"), "supply `visits`")
  # a pair whose cell types share no movable cluster is flagged, not fatal
  fit2 <- .cmp_visit_toy(n_pat = 25)
  fit2$scores$C <- fit2$scores$C[1:2, , drop = FALSE]   # C overlaps the others on one patient only
  expect_warning(ss2 <- sharing_spectrum(fit2, n_perm = 49, seed = 1, clusters = "patient", visits = "visit"),
                 "skipped")
  expect_true(all(ss2$summary$underpowered[grepl("C", ss2$summary$pair)]))
})
