## Build the comparison-vignette cache from the trained Stephenson MOFA model
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(cellprograms); library(MOFAcellulaR); library(MOFA2)
})

view_dims <- c(B_mem=848, B_naive=445, CD14_Mono=295, CD4_mem=632, CD4_naive=390, CD8_mem=528, CD8_naive=546, DC=259, NK_bright=602, NK_dim=431, Plasma=880, unconv_T=841)
f <- readRDS("/workspace/stephenson_K30_masked_fit.rds")
md <- f$sample_metadata
model <- load_model("/workspace/stephenson_mofa.hdf5")

Z <- get_factors(model)[[1]]                      # donors x factors
colnames(Z) <- paste0("MF", seq_len(ncol(Z)))
ve <- get_variance_explained(model)
r2_pf <- ve$r2_per_factor[[1]]                    # factors x views
r2_tot <- ve$r2_total[[1]]                        # per-view total R2 (named vector)

## alignment: every cellprograms program x every MOFA factor; donors per
## cell type differ (72-117), so correlate each type's programs on the
## donors that type shares with the MOFA factor matrix
r_mat <- do.call(rbind, lapply(f$cell_types, function(ct) {
  S <- f$scores[[ct]]
  d <- intersect(rownames(S), rownames(Z))
  cor(S[d, ], Z[d, ])
}))
rownames(r_mat) <- unlist(lapply(f$scores, colnames))
colnames(r_mat) <- colnames(Z)

## severity associations, both models
cmp <- readRDS("/workspace/stephenson_K30_cmp.rds")
ca <- cmp_context_association(f, cmp, "severity")
md_f <- data.frame(sample = rownames(Z),
                   severity = md$severity[match(rownames(Z), md$observation_id)])
ga <- MOFAcellulaR::get_associations(model = model, metadata = md_f,
                                     sample_id_column = "sample",
                                     test_variable = "severity",
                                     test_type = "continuous", group = FALSE)

cache <- list(
  Z = Z, r_mat = r_mat,
  r2_per_factor = r2_pf, r2_total = r2_tot,
  cmp_assoc = ca, mofa_assoc = as.data.frame(ga),
  blacklist_n = length(f$preprocessing$exclude_genes),
  n_donors = nrow(Z), n_donors_input = nrow(md),
  qc_note = "MOFAcellulaR QC chain (expression + 0.9 coverage filters) retained 118 of 119 donors; view gene counts after their expression filter: see view_dims",
  view_dims = view_dims,
  settings = list(num_factors = 15, seed = 1, spikeslab_weights = FALSE,
                  min_count = 5, min_prop = 0.25, prop_coverage = 0.9,
                  scale_factor = 1e6, prior_hvg = "masked-fit feature sets"),
  session = sessionInfo()
)
saveRDS(cache, "/workspace/repo2/vignettes/articles/mofacellular_stephenson_cache.rds")
cat("cache saved. donors:", nrow(Z), " max|r|:", round(max(abs(r_mat)), 3), "\n")
cat("top severity associations (MOFA):\n")
print(head(ga[order(ga$p.value), ], 3))
cat("top severity associations (cellprograms CMP):\n")
print(head(ca[order(ca$p), ], 3))
cat("done\n")
