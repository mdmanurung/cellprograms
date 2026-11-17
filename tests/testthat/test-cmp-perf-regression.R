# Regression guard for the find_coordinated_programs() permutation loop.
# Golden values come from the loop before the fast path (release/0.3.0 @ 0914144);
# see fixtures/make-cmp-perf-golden.R. node_p sits on a 1/(n_perm + 1) grid and
# depends on every permutation maximum, so an identical node_p pins the RNG stream
# and the maxima. Statistics may differ by floating-point noise only.

golden <- readRDS(test_path("fixtures", "cmp_perf_golden.rds"))

expect_pin_equal <- function(got, want) {
  expect_identical(got$nodes$node_p, want$nodes$node_p)
  expect_identical(got$nodes$callable, want$nodes$callable)
  expect_equal(got$nodes$node_stat, want$nodes$node_stat, tolerance = 1e-10)
  expect_identical(got$cmps[c("cmp", "cell_type", "program_id", "node_p")],
                   want$cmps[c("cmp", "cell_type", "program_id", "node_p")])
  expect_equal(got$cmps$sign, want$cmps$sign)
  expect_identical(got$edge_support[c("ct_a", "ct_b", "prog_a", "prog_b", "n", "callable", "reason")],
                   want$edge_support[c("ct_a", "ct_b", "prog_a", "prog_b", "n", "callable", "reason")])
  expect_equal(got$edge_support[c("stat", "sign")], want$edge_support[c("stat", "sign")],
               tolerance = 1e-10)
  expect_equal(got$edges, want$edges, tolerance = 1e-10)
}

cases <- .perf_cases()

for (nm in names(cases)) {
  test_that(paste0("find_coordinated_programs output is unchanged: ", nm), {
    expect_pin_equal(.perf_pin(.perf_run(cases[[nm]])), golden[[nm]])
  })
}

test_that("the golden cases exercise the paths the fast path must keep", {
  expect_gt(nrow(golden$complete$cmps), 0L)
  expect_setequal(paste(golden$complete$cmps$cell_type, golden$complete$cmps$program_id),
                  c("A A_1", "B B_2", "C C_3"))
  expect_true("insufficient_shared" %in% golden$sparse_celltype$edge_support$reason)
  expect_true("undefined_or_uncallable_residual" %in% golden$constant_program$edge_support$reason)
  expect_gt(length(unique(golden$complete$nodes$node_p)), 1L)
  expect_equal(.perf_run(cases$clusters)$settings$permutation, "cluster-within-strata-and-visit-set")
  expect_warning(do.call(find_coordinated_programs,
                         c(list(fit = cases$na_donors$fit, n_perm = 99L, seed = 1L),
                           cases$na_donors$args)), "missing donor scores")
})

test_that("the seed decides the result and the stream is stable across calls", {
  a <- .perf_pin(.perf_run(cases$complete))
  expect_identical(.perf_pin(.perf_run(cases$complete)), a)
  other <- .perf_pin(.perf_run(cases$complete, seed = 2L))
  expect_false(identical(other$nodes$node_p, a$nodes$node_p))
})

test_that("assess_cmp_stability output is unchanged", {
  got <- .perf_stability_pin()
  want <- golden$stability
  expect_equal(got$summary, want$summary, tolerance = 1e-10)
  expect_equal(got$members, want$members, tolerance = 1e-10)
  expect_equal(got$replicates, want$replicates, tolerance = 1e-10)
})

# ---- internals the fast path has to match -------------------------------------------------

test_that(".cmp_perm draws for every cell type, in names(Z) order, even those in no pair", {
  ids <- sprintf("d%02d", 1:30)
  mk <- function(rows) matrix(rnorm(length(rows) * 2), length(rows), 2,
                              dimnames = list(rows, c("x1", "x2")))
  set.seed(1)
  Z <- list(A = mk(ids), B = mk(ids), C = mk(sprintf("e%02d", 1:10)))
  blocks <- stats::setNames(rep("s", 40), c(ids, rownames(Z$C)))
  set.seed(5); .cmp_perm(Z, blocks); with_c <- runif(1)
  set.seed(5); .cmp_perm(Z[c("A", "B")], blocks); without_c <- runif(1)
  expect_false(isTRUE(all.equal(with_c, without_c)))
})

test_that(".cmp_edges keeps a finite stat when a program has no support row (callable NA)", {
  set.seed(3)
  ids <- sprintf("d%02d", 1:30)
  Za <- matrix(rnorm(60), 30, 2, dimnames = list(ids, c("A_1", "A_2")))
  Zb <- matrix(rnorm(90), 30, 3, dimnames = list(ids, c("B_1", "B_2", "B_3")))
  attr(Za, "program_support") <- data.frame(program_id = c("A_1", "A_2"), callable = c(TRUE, TRUE))
  attr(Zb, "program_support") <- data.frame(program_id = c("B_1", "B_2"), callable = c(TRUE, TRUE))
  e <- .cmp_edges(list(A = Za, B = Zb), list(list(a = "A", b = "B", ids = ids, min_shared = 4L)))
  b3 <- e[e$prog_b == "B_3", ]
  expect_true(all(is.na(b3$callable)))
  expect_true(all(is.finite(b3$stat)))
})

test_that(".cmp_edges clamps |r| at 0.999999 and masks n <= 3", {
  ids <- sprintf("d%02d", 1:30)
  a <- matrix(seq_len(30) + rnorm(30, sd = 0.01), 30, 1, dimnames = list(ids, "A_1"))
  b <- cbind(B_1 = a[, 1], B_2 = -a[, 1])
  pr <- list(a = "A", b = "B", ids = ids, min_shared = 4L)
  e <- .cmp_edges(list(A = a, B = b), list(pr))
  expect_equal(e$stat, rep(atanh(0.999999) * sqrt(27), 2), tolerance = 1e-6)
  expect_equal(e$sign, c(1, -1))
  few <- .cmp_edges(list(A = a, B = b), list(modifyList(pr, list(ids = ids[1:3], min_shared = 1L))))
  expect_true(all(is.na(few$stat)))
  expect_false(any(few$callable))
})

test_that(".cmp_node_stat ignores non-finite edge stats and sums over cell-type partners", {
  edges <- data.frame(ct_a = "A", ct_b = "B", prog_a = c("A_1", "A_2"), prog_b = "B_1")
  idx <- .cmp_node_index(edges)
  expect_equal(unname(.cmp_node_stat(c(NA, 5), idx, 2)), c(0, 9, 9))
  expect_equal(unname(.cmp_node_stat(c(Inf, 1), idx, 2)), c(0, 0, 0))
  expect_equal(unname(.cmp_node_stat(c(1, 1), idx, -1)), c(4, 4, 4))
})
