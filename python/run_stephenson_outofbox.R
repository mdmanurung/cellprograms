# Out-of-box R package pipeline on Stephenson log-norm+HVG pseudobulk.
# Pure defaults: point_laplace, max_factors=10, var_type=1, backfit, seed=1.
# Emits the adapter-equivalent representation (z-scaled score blocks, concat).
args <- commandArgs(trailingOnly = TRUE)
arm_dir <- args[1]; out_csv <- args[2]
suppressMessages(library(cellprograms))
suppressMessages(library(matrixStats))

files <- sort(list.files(arm_dir, pattern = "^lognorm_hvg__.*\\.csv$", full.names = TRUE))
pb <- lapply(files, function(f) as.matrix(read.csv(f, row.names = 1, check.names = FALSE)))
names(pb) <- sub("^lognorm_hvg__|\\.csv$", "", basename(files))
cat("cell types:", length(pb), "\n", sep = " ")

cpd <- as_cell_program_data(pb)
t0 <- Sys.time()
fit <- canonicalize_programs(fit_celltype_programs(cpd))  # pure defaults, seed=1
cat("fit: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n", sep = "")

blocks <- lapply(names(fit$scores), function(ct) {
  S <- fit$scores[[ct]]
  if (is.null(dim(S)) || ncol(S) == 0) return(NULL)
  S <- scale(as.matrix(S))  # adapter block_scale="z" equivalent
  colnames(S) <- paste0(ct, "__", seq_len(ncol(S)))
  S
})
names(blocks) <- names(fit$scores)
blocks <- blocks[!vapply(blocks, is.null, logical(1))]
all_samp <- sort(unique(unlist(lapply(blocks, rownames))))
rep <- do.call(cbind, lapply(blocks, function(B) {
  M <- matrix(0, nrow = length(all_samp), ncol = ncol(B),
              dimnames = list(all_samp, colnames(B)))
  M[rownames(B), ] <- B
  M
}))
rep[is.na(rep)] <- 0
write.csv(rep, out_csv)
cat("representation:", nrow(rep), "x", ncol(rep), "\n", sep = "")
