## run_mofa_stephenson_excl.R -- MOFAcellulaR sensitivity run: expanded blacklist.
## Same chain as run_mofa_stephenson.R but (a) gene panels from the
## exclude_artifacts=TRUE refit and (b) the full 4-set artifact blacklist vetoed
## via the phantom-origin view. Train-only; the comparison cache is built
## separately from the saved model.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(cellprograms)
  library(MOFAcellulaR)
  library(MOFA2)
})

f <- readRDS("/workspace/stephenson_K30_masked_fit_excl.rds")
cts <- f$cell_types
shared <- f$shared_genes
md <- f$sample_metadata

counts <- do.call(cbind, lapply(cts, function(ct) {
  t(f$pseudobulk[[ct]][, shared, drop = FALSE])
}))
ct_of <- rep(cts, times = vapply(f$pseudobulk, nrow, 1L))
donor_of <- unlist(lapply(f$pseudobulk, rownames), use.names = FALSE)
colnames(counts) <- paste(ct_of, donor_of, sep = ".")

coldata <- do.call(rbind, lapply(cts, function(ct) {
  d <- md[match(rownames(f$pseudobulk[[ct]]), md$observation_id), , drop = FALSE]
  d$cell_type <- ct
  d$cell_counts <- 0L
  d
}))
rownames(coldata) <- colnames(counts)

## expanded blacklist (all 4 artifact sets), restricted to the view universe
blacklist <- artifact_genes(shared)
cat("shared genes:", length(shared), "| vetoed:", length(blacklist), "\n")
flush.console()

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
flush.console()

multiview_dat <- MOFAcellulaR::pb_dat2MOFA(pb_dat_list = ct_list_se,
                                           sample_column = "observation_id")
multiview_dat <- MOFA2::create_mofa(multiview_dat)

mo_opts <- MOFA2::get_default_model_options(multiview_dat)
mo_opts$num_factors <- 15
mo_opts$spikeslab_weights <- FALSE
tr_opts <- MOFA2::get_default_training_options(multiview_dat)
tr_opts$seed <- 1
multiview_dat <- MOFA2::prepare_mofa(multiview_dat,
                                     data_options = MOFA2::get_default_data_options(multiview_dat),
                                     model_options = mo_opts, training_options = tr_opts)
cat("training MOFA2 (15 factors, 12 views)...\n"); flush.console()
model <- MOFA2::run_mofa(multiview_dat, outfile = "/workspace/stephenson_mofa_excl.hdf5")
cat("MODEL SAVED\n")
