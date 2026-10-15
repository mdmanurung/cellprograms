## 05_build_cache.R -- build the comparison-vignette cache from a trained MOFA
## model + the matching cellprograms fit and CMP result. Produces the small
## object the vignette reads (donor-score alignment, severity associations,
## per-view R2, view dims, settings).
##
## Usage: Rscript 05_build_cache.R <fit_rds> <cmp_rds> <mofa_hdf5> <out_rds>
##   e.g. Rscript 05_build_cache.R results/stephenson_K30_masked_fit.rds \
##          results/stephenson_K30_cmp.rds results/stephenson_mofa.hdf5 \
##          vignettes/articles/mofacellular_stephenson_cache.rds
##
## view_dims (genes per view after the MOFAcellulaR QC chain) are read from the
## 03_mofa.R log; pass them via the VIEW_DIMS env var as a comma-separated
## string in cell-type order, or edit the default below.

suppressPackageStartupMessages({
  library(cellprograms); library(MOFAcellulaR); library(MOFA2)
})

a <- commandArgs(TRUE)
stopifnot("usage: 05_build_cache.R <fit_rds> <cmp_rds> <mofa_hdf5> <out_rds>" = length(a) >= 4)
fit_rds <- a[1]; cmp_rds <- a[2]; mofa_hdf5 <- a[3]; out_rds <- a[4]

f <- readRDS(fit_rds)
md <- f$sample_metadata
model <- load_model(mofa_hdf5)

vd_env <- Sys.getenv("VIEW_DIMS", "")
view_dims <- if (nzchar(vd_env)) {
  as.integer(strsplit(vd_env, ",")[[1]])
} else {
  c(B_mem = 848, B_naive = 445, CD14_Mono = 295, CD4_mem = 632, CD4_naive = 390,
    CD8_mem = 528, CD8_naive = 546, DC = 259, NK_bright = 602, NK_dim = 431,
    Plasma = 880, unconv_T = 841)   # default-setting run; override via VIEW_DIMS
}

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

cmp <- readRDS(cmp_rds)
ca <- cmp_context_association(f, cmp, "severity")
md_f <- data.frame(sample = rownames(Z),
                   severity = md$severity[match(rownames(Z), md$observation_id)])
ga <- MOFAcellulaR::get_associations(model = model, metadata = md_f,
                                     sample_id_column = "sample",
                                     test_variable = "severity",
                                     test_type = "continuous", group = FALSE)

is_excl <- isTRUE(f$preprocessing$exclude_artifacts)
cache <- list(
  Z = Z, r_mat = r_mat, r2_per_factor = r2_pf, r2_total = r2_tot,
  cmp_assoc = ca, mofa_assoc = as.data.frame(ga),
  blacklist_n = length(f$preprocessing$exclude_genes),
  n_donors = nrow(Z), n_donors_input = nrow(md),
  qc_note = sprintf(paste0("MOFAcellulaR QC chain (expression + 0.9 coverage ",
                           "filters) retained %d of %d donors; view gene ",
                           "counts after their expression filter: see view_dims"),
                    nrow(Z), nrow(md)),
  view_dims = view_dims,
  settings = list(num_factors = 15, seed = 1, spikeslab_weights = FALSE,
                  min_count = 5, min_prop = 0.25, prop_coverage = 0.9,
                  scale_factor = 1e6,
                  prior_hvg = if (is_excl) "exclude_artifacts refit feature sets" else
                    "masked-fit feature sets"),
  session = sessionInfo()
)
if (is_excl) cache$blacklist_n_default <- 274L   # default-arm blacklist size
saveRDS(cache, out_rds)
cat("cache saved:", out_rds, "| donors:", nrow(Z),
    "| max|r|:", round(max(abs(r_mat)), 3), "\n")
