## build_cache_excl.R -- comparison cache for the exclude_artifacts=TRUE setting.
## Same structure as the default-setting cache; adds the old-vs-new blacklist
## sizes so the vignette sensitivity subsection can state the delta exactly.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(cellprograms); library(MOFAcellulaR); library(MOFA2)
})

view_dims <- c(B_mem = 846, B_naive = 433, CD14_Mono = 286, CD4_mem = 621,
               CD4_naive = 376, CD8_mem = 515, CD8_naive = 534, DC = 250,
               NK_bright = 595, NK_dim = 419, Plasma = 874, unconv_T = 822)

f <- readRDS("/workspace/stephenson_K30_masked_fit_excl.rds")
md <- f$sample_metadata
model <- load_model("/workspace/stephenson_mofa_excl.hdf5")

Z <- get_factors(model)[[1]]
colnames(Z) <- paste0("MF", seq_len(ncol(Z)))
ve <- get_variance_explained(model)
r2_pf <- ve$r2_per_factor[[1]]
r2_tot <- ve$r2_total[[1]]

## alignment: per cell type on the donors that type shares with Z
r_mat <- do.call(rbind, lapply(f$cell_types, function(ct) {
  S <- f$scores[[ct]]
  d <- intersect(rownames(S), rownames(Z))
  cor(S[d, ], Z[d, ])
}))
rownames(r_mat) <- unlist(lapply(f$scores, colnames))
colnames(r_mat) <- colnames(Z)

cmp <- readRDS("/workspace/stephenson_K30_cmp_excl.rds")
ca <- cmp_context_association(f, cmp, "severity")
md_f <- data.frame(sample = rownames(Z),
                   severity = md$severity[match(rownames(Z), md$observation_id)])
ga <- MOFAcellulaR::get_associations(model = model, metadata = md_f,
                                     sample_id_column = "sample",
                                     test_variable = "severity",
                                     test_type = "continuous", group = FALSE)

all_genes <- unique(unlist(lapply(f$pseudobulk, colnames)))
cache <- list(
  Z = Z, r_mat = r_mat,
  r2_per_factor = r2_pf, r2_total = r2_tot,
  cmp_assoc = ca, mofa_assoc = as.data.frame(ga),
  blacklist_n = length(f$preprocessing$exclude_genes),
  blacklist_n_default = 274L,
  n_donors = nrow(Z), n_donors_input = nrow(md),
  view_dims = view_dims,
  settings = list(num_factors = 15, seed = 1, spikeslab_weights = FALSE,
                  min_count = 5, min_prop = 0.25, prop_coverage = 0.9,
                  scale_factor = 1e6, prior_hvg = "exclude_artifacts refit feature sets",
                  exclude_artifacts = TRUE),
  session = sessionInfo()
)
saveRDS(cache, "/workspace/repo2/vignettes/articles/mofacellular_stephenson_cache_excl.rds")
cat("cache saved. donors:", nrow(Z), " max|r|:", round(max(abs(r_mat)), 3), "\n")
cat("top severity (MOFA):\n"); print(head(ga[order(ga$p.value), ], 3))
cat("top severity (CMP):\n"); print(head(ca[order(ca$p), ], 3))
cat("done\n")
