# MOFAcellulaR fit via standalone Rscript (avoids rpy2/reticulate
# Python-in-R-in-Python conflict with mofapy2).
# Usage: Rscript mofa_fit.R <mat.rds> <coldata.rds> <tag> <n_factors> <seed> <hvg:0|1> <out.csv>
suppressMessages({
  .libPaths(c("/workspace/.Rlib", .libPaths()))
  library(MOFAcellulaR); library(MOFA2); library(scran)
})
args <- commandArgs(trailingOnly = TRUE)
mat <- readRDS(args[1]); coldata <- readRDS(args[2])
tag <- args[3]; n_factors <- as.integer(args[4]); seed <- as.integer(args[5])
hvg <- as.integer(args[6]); outfile <- args[7]

assays <- if (tag == "A") list(counts = mat, logcounts = mat) else list(counts = mat)
se <- SummarizedExperiment::SummarizedExperiment(
  assays = assays, colData = S4Vectors::DataFrame(coldata))
pb_list <- filt_profiles(se, cts = NULL, ncells = 50,
                         counts_col = "cell_counts", ct_col = "cell_type")
if (tag != "A") {
  pb_list <- tmm_trns(pb_list, scale_factor = 1e6)
  if (hvg == 1) {
    pb_list <- lapply(pb_list, function(x) {
      keep <- scran::getTopHVGs(scran::modelGeneVar(assay(x, "logcounts")), n = 2000)
      x[keep, ]
    })
  }
}
mofa_list <- pb_dat2MOFA(pb_list, sample_column = "donor_id")
mofa <- create_mofa(mofa_list)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- n_factors
train_opts <- get_default_training_options(mofa)
train_opts$seed <- seed
train_opts$convergence_mode <- "fast"
mofa <- prepare_mofa(mofa, data_options = get_default_data_options(mofa),
                     model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa, outfile = sprintf("/workspace/combat_benchmark/mofa_%s.hdf5", tag))
fac <- get_factors(mofa, factors = "all")[[1]]
out <- data.frame(sample = rownames(fac), as.matrix(fac), check.names = FALSE)
write.csv(out, outfile, row.names = FALSE)
cat("MOFA fit done:", tag, nrow(fac), "samples x", ncol(fac), "factors\n")
