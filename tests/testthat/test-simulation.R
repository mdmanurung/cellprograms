test_that("mixed scenario shares activity within the labelled shared programs", {
  s <- simulate_pseudobulk("mixed", n_samples = 30, n_genes = 200, program_size = 10, seed = 2)
  z_of <- function(label) {
    out <- list()
    for (ct in names(s$truth)) for (p in s$truth[[ct]]) if (p$label == label) out[[ct]] <- p$z
    out
  }
  ifn <- z_of("global_IFN")
  expect_setequal(names(ifn), c("B", "CD4", "NK", "Mono"))
  for (z in ifn[-1]) expect_equal(z, ifn[[1]])
  cn <- z_of("CD8NK_shared")
  expect_equal(cn$CD8, cn$NK)
})
