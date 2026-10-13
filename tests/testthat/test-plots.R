test_that("plot functions return ggplot objects on a standard fit", {
  skip_if_not_installed("ggplot2")
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
