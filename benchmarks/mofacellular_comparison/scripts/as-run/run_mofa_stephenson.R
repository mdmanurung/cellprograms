## MOFAcellulaR on the Stephenson pseudobulk, same gene panels as the masked
## cellprograms fits. Produces the small cache used by the comparison vignette.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(cellprograms)
  library(MOFAcellulaR)
  library(MOFA2)
})

f <- readRDS("/workspace/stephenson_K30_masked_fit.rds")
cts <- f$cell_types
shared <- f$shared_genes
md <- f$sample_metadata

## -- build the (genes x pseudobulk profiles) matrix + colData the way
##    create_init_exp expects: one SE, coldata keyed by profile
counts <- do.call(cbind, lapply(cts, function(ct) {
  t(f$pseudobulk[[ct]][, shared, drop = FALSE])          # genes x donors
}))
ct_of <- rep(cts, times = vapply(f$pseudobulk, nrow, 1L))
donor_of <- unlist(lapply(f$pseudobulk, rownames), use.names = FALSE)
colnames(counts) <- paste(ct_of, donor_of, sep = ".")

coldata <- do.call(rbind, lapply(cts, function(ct) {
  d <- md[match(rownames(f$pseudobulk[[ct]]), md$observation_id), , drop = FALSE]
  d$cell_type <- ct
  d$cell_counts <- 0L                                    # not available; no filter
  d
}))
rownames(coldata) <- colnames(counts)
stopifnot(all(coldata$observation_id %in% md$observation_id))

## -- same blacklist as the masked fits (fit$preprocessing$exclude_genes),
##    removed from every view via a phantom origin view
blacklist <- f$preprocessing$exclude_genes
cat("blacklist genes:", length(blacklist), "\n")

cat("running MOFAcellulaR QC chain...\n")
ct_list_se <- MOFAcellulaR::create_init_exp(counts = counts, coldata = coldata) %>%
  MOFAcellulaR::filt_profiles(pb_dat = ., cts = cts, ncells = 0,
                              counts_col = "cell_counts", ct_col = "cell_type") %>%
  MOFAcellulaR::filt_views_bysamples(pb_dat_list = ., nsamples = 2) %>%
  MOFAcellulaR::filt_gex_byexpr(pb_dat_list = ., min.count = 5, min.prop = 0.25) %>%
  MOFAcellulaR::filt_views_bygenes(pb_dat_list = ., ngenes = 15) %>%
  MOFAcellulaR::filt_samples_bycov(pb_dat_list = ., prop_coverage = 0.9) %>%
  MOFAcellulaR::tmm_trns(pb_dat_list = ., scale_factor = 1e6) %>%
  MOFAcellulaR::filt_gex_byhvg(pb_dat_list = .,
                               prior_hvg = f$preprocessing$feature_sets,
                               var.threshold = 0) %>%
  MOFAcellulaR::filt_gex_bybckgrnd(pb_dat_list = .,
                                   prior_mrks = list(phantom_origin = blacklist)) %>%
  MOFAcellulaR::filt_views_bygenes(pb_dat_list = ., ngenes = 15)

cat("views after QC (genes per view):\n")
view_dims <- sapply(ct_list_se, nrow)
print(view_dims)

multiview_dat <- MOFAcellulaR::pb_dat2MOFA(pb_dat_list = ct_list_se,
                                           sample_column = "observation_id")
multiview_dat <- MOFA2::create_mofa(multiview_dat)

## -- train
mo_opts <- MOFA2::get_default_model_options(multiview_dat)
mo_opts$num_factors <- 15
mo_opts$spikeslab_weights <- FALSE
tr_opts <- MOFA2::get_default_training_options(multiview_dat)
tr_opts$seed <- 1
multiview_dat <- MOFA2::prepare_mofa(multiview_dat,
                                     data_options = MOFA2::get_default_data_options(multiview_dat),
                                     model_options = mo_opts, training_options = tr_opts)
cat("training MOFA2 (15 factors, 12 views)...\n")
model <- MOFA2::run_mofa(multiview_dat, outfile = "/workspace/stephenson_mofa.hdf5")

## -- cache the small objects the vignette needs
Z <- MOFA2::get_factors(model)$group1
colnames(Z) <- paste0("MF", seq_len(ncol(Z)))
r2_list <- MOFA2::get_variance_explained(model)$R2
r2 <- if (is.list(r2_list)) r2_list[[1]] else r2_list     # group1: views x factors
if (nrow(r2) != length(view_dims)) r2 <- t(r2)            # force rows = views
r2_df <- as.data.frame(as.table(r2))
colnames(r2_df) <- c("view", "factor", "r2")

## alignment: every cellprograms program x every MOFA factor
progs <- do.call(cbind, f$scores)
common <- intersect(rownames(progs), rownames(Z))
r_mat <- cor(progs[common, ], Z[common, ])

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
  Z = Z, r_mat = r_mat, r2_df = r2_df, view_dims = view_dims,
  cmp_assoc = ca, mofa_assoc = as.data.frame(ga),
  blacklist_n = length(blacklist), n_donors = nrow(Z),
  settings = list(num_factors = 15, seed = 1, spikeslab_weights = FALSE,
                  min_count = 5, min_prop = 0.25, prop_coverage = 0.9,
                  scale_factor = 1e6, prior_hvg = "masked-fit feature sets"),
  session = sessionInfo()
)
saveRDS(cache, "/workspace/repo2/vignettes/articles/mofacellular_stephenson_cache.rds")
cat("cache saved. dim(Z):", dim(Z), " max|r|:", max(abs(r_mat)), "\n")
cat("done\n")
