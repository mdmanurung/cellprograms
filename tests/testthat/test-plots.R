test_that("plot functions return ggplot objects on a standard fit", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("mashr") # program_sharing() below needs mashr (Suggests)
  sim <- simulate_pseudobulk(n_samples = 30, n_genes = 300, program_size = 30, seed = 1)
  z <- sim$truth$B[[1]]$z > median(sim$truth$B[[1]]$z)
  md <- data.frame(observation_id = paste0("obs", seq_along(z)),
                   severity = factor(z, c(FALSE, TRUE), c("low", "high")))
  x <- as_cell_program_data(sim$pseudobulk, sample_metadata = md)
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 5, seed = 1))
  sh <- program_sharing(fit, n_perm = 200, seed = 1)
  sp <- sharing_spectrum(fit, n_perm = 200)
  cm <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
  fit_st <- assess_program_stability(fit, n_boot = 5, seed = 1)

  is_gg <- function(p) inherits(p, "ggplot")
  expect_true(is_gg(plot_program_scores(fit)))
  expect_true(is_gg(plot_program_scores(fit, annotation = "severity")))
  expect_true(is_gg(plot_top_genes(fit, c("B_1", "CD4_1"), n = 8)))
  expect_true(is_gg(plot_variance_explained(fit)))
  expect_true(is_gg(plot_sharing(sh)))
  expect_true(is_gg(plot_sharing_spectrum(sp)))
  expect_true(is_gg(plot_artifact_flags(flag_artifact_genes(fit))))
  expect_true(is_gg(plot_cmp_network(cm)))
  expect_true(is_gg(plot_cmp_nodes(cm)))
  expect_true(is_gg(plot_context_association(cmp_context_association(fit, cm, "severity"))))
  expect_true(is_gg(plot_alignment(principal_angles(fit, "B", "CD4", n_perm = 200))))
  expect_true(is_gg(plot_alignment(align_programs(fit, n_perm = 200))))
  expect_true(is_gg(plot_stability(fit_st)))
})

test_that("plot functions degrade gracefully on degenerate inputs", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_pseudobulk(n_samples = 30, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 5, seed = 1))

  ## K = 0 cell type
  fit_k0 <- fit
  ct <- "CD8"
  fit_k0$scores[[ct]] <- fit_k0$scores[[ct]][, 0, drop = FALSE]
  fit_k0$loadings[[ct]] <- fit_k0$loadings[[ct]][, 0, drop = FALSE]
  expect_true(inherits(plot_program_scores(fit_k0), "ggplot"))
  expect_true(inherits(plot_variance_explained(fit_k0), "ggplot"))
  expect_true(inherits(plot_top_genes(fit_k0, "CD8_1"), "ggplot"))

  ## all cell types K = 0
  fit_all0 <- fit
  for (c0 in names(fit_all0$scores)) {
    fit_all0$scores[[c0]] <- fit_all0$scores[[c0]][, 0, drop = FALSE]
    fit_all0$loadings[[c0]] <- fit_all0$loadings[[c0]][, 0, drop = FALSE]
  }
  expect_message(p0 <- plot_program_scores(fit_all0), "No programs")
  expect_true(inherits(p0, "ggplot"))

  ## zero-edge CMP result
  cm0 <- find_coordinated_programs(fit, n_perm = 200, seed = 1, link_threshold = 1e9)
  expect_message(pn <- plot_cmp_network(cm0), "No significant coordination edges")
  expect_true(inherits(pn, "ggplot"))
  expect_true(inherits(plot_cmp_nodes(cm0), "ggplot"))

  ## empty inputs
  expect_message(pe <- plot_sharing(list(matches = data.frame(), summary = data.frame()),
                                    q_alpha = 0.05), "No cross-cell-type")
  expect_true(inherits(pe, "ggplot"))
  expect_message(ps <- plot_stability(fit), "assess_program_stability")
  expect_true(inherits(ps, "ggplot"))
  expect_message(pa <- plot_context_association(data.frame()), "No context associations")
  expect_true(inherits(pa, "ggplot"))
})

test_that("plot_cmp_network falls back without igraph", {
  skip_if_not_installed("ggplot2")
  skip_if(!requireNamespace("igraph", quietly = TRUE))
  sim <- simulate_pseudobulk(n_samples = 30, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 5, seed = 1))
  cm <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
  p1 <- plot_cmp_network(cm, seed = 1)
  expect_true(inherits(p1, "ggplot"))
  ## max_edges truncation path
  p2 <- plot_cmp_network(cm, max_edges = 1)
  expect_true(inherits(p2, "ggplot"))
})

test_that("plot_cmp_membership classifies, canonicalizes labels, and degrades", {
  sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
  x <- as_cell_program_data(sim$pseudobulk)
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
  cmp_d <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
  cmp_f <- find_coordinated_programs(fit, n_perm = 200, seed = 1, ambient = "flag")
  p <- plot_cmp_membership(list(drop = cmp_d, flag = cmp_f))
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  ## degenerate: nothing called -> still a valid plot
  cmp0 <- find_coordinated_programs(fit, n_perm = 200, seed = 2, tau = 1e5)
  p0 <- plot_cmp_membership(list(none = cmp0))
  expect_s3_class(p0, "ggplot")
  expect_silent(ggplot2::ggplot_build(p0))
  ## input validation
  expect_error(plot_cmp_membership(list(cmp_d)), "named list")
  expect_error(plot_cmp_membership(list(a = fit)), "cmp_result")
})

test_that("plot_cmp_network handles results computed with ambient = 'none' (no ambient_score column)", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_pseudobulk(n_samples = 30, n_genes = 300, program_size = 30, seed = 1)
  md <- data.frame(observation_id = paste0("obs", seq_len(30)),
                   site = rep(c("s1", "s2"), length.out = 30))
  x <- as_cell_program_data(sim$pseudobulk, sample_metadata = md)
  fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 5, seed = 1))
  cm_none <- find_coordinated_programs(fit, strata = "site", n_perm = 99,
                                       ambient = "none", seed = 1)
  expect_null(cm_none$nodes$ambient_score)   # the column is genuinely absent
  p <- plot_cmp_network(cm_none)             # regressed: setNames length mismatch
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_cmp_network marks ambient-flagged programs with the hollow-triangle shape", {
  skip_if_not_installed("ggplot2")
  ids <- c("A_1", "B_1", "C_1")
  cmp <- structure(list(
    edges = data.frame(prog_a = c("A_1", "A_1"), prog_b = c("B_1", "C_1"),
                       stat = c(4, 3), sign = c(1, -1), stringsAsFactors = FALSE),
    nodes = data.frame(cell_type = c("A", "B", "C"), program_id = ids,
                       node_stat = c(5, 4, 3), ambient_score = c(0.9, 0.1, 0.1),
                       stringsAsFactors = FALSE),
    cmps = data.frame(cmp = 1L, cell_type = c("A", "B", "C"), program_id = ids,
                      stringsAsFactors = FALSE),
    settings = list(ambient_threshold = 0.3)), class = "cmp_result")
  p <- plot_cmp_network(cmp, seed = 1)
  pts <- ggplot2::ggplot_build(p)$data[[2]]
  ## shape 24 = ambient; exactly one of the three nodes (A_1) is flagged
  expect_equal(sum(pts$shape == 24), 1L)
})

test_that("plot_alignment reports significance from align_programs() output", {
  skip_if_not_installed("ggplot2")
  mk <- function(p) structure(list(
    program_weights = list(a = matrix(c(1, 0), 2, 1, dimnames = list(c("A_1", "A_2"), "d1")),
                           b = matrix(c(0, 1), 2, 1, dimnames = list(c("B_1", "B_2"), "d1"))),
    p_pair = p, underpowered = FALSE), class = "principal_angles")
  x <- structure(list(B = mk(0.001), C = mk(0.4)), reference = "A")
  p <- plot_alignment(x)
  expect_match(p$labels$subtitle, "1/2 pairs significant")
  expect_match(p$labels$subtitle, "unadjusted p")
})
