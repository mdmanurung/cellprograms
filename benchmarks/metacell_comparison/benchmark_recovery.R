# Stage 4: downstream program recovery from metacell counts.
# Usage: Rscript benchmark_recovery.R <fixture_dir> <mem_dir> <scale> <seed> <method> <N> <out_csv>
# Seam (same for both methods): metacell counts -> dominant (cell type, donor) per
# metacell -> pseudobulk counts per (type, donor) -> as_cell_program_data ->
# fit_celltype_programs -> detection (loading-vs-planted-signature correlation)
# + fidelity (metacell program scores vs planted donor activity).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(cellprograms)
})

args <- commandArgs(trailingOnly = TRUE)
fx_dir <- args[1]; mem_dir <- args[2]; scale <- args[3]; seed <- as.integer(args[4])
method <- args[5]; N <- as.integer(args[6]); out_csv <- args[7]
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(fx_dir, paste0("counts_", scale, "_s", seed, ".rds")))
meta <- readRDS(file.path(fx_dir, paste0("meta_", scale, "_s", seed, ".rds")))
program_genes <- readRDS(file.path(fx_dir, paste0("programs_", scale, "_s", seed, ".rds")))
donor_activity <- readRDS(file.path(fx_dir, paste0("activity_", scale, "_s", seed, ".rds")))
mem <- read.csv(file.path(mem_dir, sprintf("membership_%s_%s_s%d_N%d.csv",
                                           method, scale, seed, N)),
                stringsAsFactors = FALSE)
mem_vec <- mem$membership; names(mem_vec) <- mem$cell_id

# --- aggregate counts per metacell (sum, matching SEACells default mode) ---
uid <- unique(mem_vec)
ind <- sparseMatrix(i = match(mem_vec, uid), j = seq_along(mem_vec),
                    dims = c(length(uid), length(mem_vec)))
mc_counts <- as.matrix(counts %*% t(ind))   # genes x metacells
colnames(mc_counts) <- as.character(uid)
rownames(mc_counts) <- rownames(counts)

# --- dominant cell type + donor per metacell ---
dom <- function(labels) {
  tab <- table(labels)
  names(tab)[which.max(tab)]
}
mc_type <- sapply(uid, function(mc) dom(meta$cell_type[match(names(mem_vec)[mem_vec == mc],
                                                             meta$cell_id)]))
mc_donor <- sapply(uid, function(mc) dom(meta$donor[match(names(mem_vec)[mem_vec == mc],
                                                          meta$cell_id)]))

# --- pseudobulk per (cell type, donor) ---
grp <- paste(mc_type, mc_donor, sep = "_")
pb <- lapply(split(colnames(mc_counts), grp), function(mcs) {
  as.matrix(rowSums(mc_counts[, mcs, drop = FALSE]))
})
all_g <- rownames(mc_counts)
pb_mat <- do.call(cbind, pb)
colnames(pb_mat) <- names(pb)

# observation metadata: donor per pseudobulk column
obs <- data.frame(
  cell_type = sub("_D.*$", "", colnames(pb_mat)),
  donor = sub("^.*_", "", colnames(pb_mat)),
  stringsAsFactors = FALSE, row.names = colnames(pb_mat)
)

# per-type pseudobulk: rows = donor observation IDs, columns = genes
pb_list <- lapply(sort(unique(obs$cell_type)), function(ct) {
  cols <- rownames(obs)[obs$cell_type == ct]
  m <- t(pb_mat[, cols, drop = FALSE])   # donors x genes
  rownames(m) <- obs$donor[match(cols, colnames(pb_mat))]
  m
})
names(pb_list) <- sort(unique(obs$cell_type))
x <- as_cell_program_data(pseudobulk = pb_list)

# max_factors must respect residual df: n_donors - 1 per type (min across types).
# At coarse N the dominant-donor design is thin (few distinct donors per type),
# so fewer factors are identifiable — part of the granularity signal.
min_donors <- min(vapply(pb_list, nrow, integer(1)))
max_factors_use <- max(1L, min(6L, min_donors - 1L))

# Downstream cost: wall time and peak allocation of the fit itself.
# gc(reset=TRUE) zeroes the per-generation maxima; the following gc() reports
# the peak (Max Used) reached during the fit — a proxy for peak memory.
gc(reset = TRUE, verbose = FALSE)
fit_t_start <- Sys.time()
fit <- tryCatch({
  f <- fit_celltype_programs(x, loading_prior = "point_laplace",
                             max_factors = max_factors_use,
                             features = "variable", n_variable_genes = 1000,
                             feature_scale = "simulated", seed = 1)
  canonicalize_programs(f)
}, error = function(e) {
  cat(sprintf("RECOVERY FIT FAIL %s %s s%d N=%d: %s\n", method, scale, seed, N,
              conditionMessage(e)))
  NULL
})
fit_time_min <- as.numeric(difftime(Sys.time(), fit_t_start, units = "mins"))
gc_res <- gc(verbose = FALSE)
# gc() returns duplicated "(Mb)" colnames; column 6 is "max used (Mb)"
fit_mem_gb <- sum(gc_res[, 6]) / 1024
if (is.null(fit)) {
  write.csv(data.frame(method = method, scale = scale, seed = seed, target_N = N,
                       n_metacells = length(uid), fit_failed = TRUE,
                       fit_time_min = fit_time_min, fit_mem_gb = fit_mem_gb),
            out_csv, row.names = FALSE)
  quit(save = "no", status = 0)
}

# --- detection: max Spearman(recovered loading vector, planted signature) per type ---
# program_loadings: list per cell type of genes x programs
loads <- program_loadings(fit)
detect <- vapply(names(program_genes), function(ct) {
  if (!ct %in% names(loads)) return(NA_real_)
  L <- loads[[ct]]
  if (is.null(dim(L)) || ncol(L) == 0) return(NA_real_)
  sig <- as.integer(rownames(L) %in% program_genes[[ct]])   # fitted genes only
  cors <- apply(L, 2, function(l) {
    r <- suppressWarnings(cor(l, sig, method = "spearman"))
    if (is.na(r)) 0 else r
  })
  max(cors)
}, numeric(1))

# --- fidelity: metacell-level program scores vs planted donor activity ---
# program_scores: per cell type, donors x programs; planted activity per donor
scores <- program_scores(fit, format = "list")
fidel <- vapply(names(program_genes), function(ct) {
  if (!ct %in% names(scores)) return(NA_real_)
  S <- scores[[ct]]
  act <- donor_activity[[ct]]
  donors_ct <- intersect(rownames(S), names(act))
  if (length(donors_ct) < 5 || ncol(S) == 0) return(NA_real_)
  cors <- apply(S[donors_ct, , drop = FALSE], 2, function(s) {
    r <- suppressWarnings(cor(s, act[donors_ct], method = "spearman"))
    if (is.na(r)) 0 else abs(r)
  })
  max(cors)
}, numeric(1))

out <- data.frame(
  method = method, scale = scale, seed = seed, target_N = N,
  n_metacells = length(uid), fit_failed = FALSE,
  fit_time_min = fit_time_min, fit_mem_gb = fit_mem_gb,
  detection_T1 = detect[["T1"]], detection_T2 = detect[["T2"]],
  detection_T3 = detect[["T3"]], detection_T4 = detect[["T4"]],
  detection_T5 = detect[["T5"]],
  detection_mean = mean(detect, na.rm = TRUE),
  fidelity_T1 = fidel[["T1"]], fidelity_T2 = fidel[["T2"]],
  fidelity_T3 = fidel[["T3"]], fidelity_T4 = fidel[["T4"]],
  fidelity_T5 = fidel[["T5"]],
  fidelity_mean = mean(fidel, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("RECOVERY OK %s %s s%d N=%d: detect=%.3f fid=%.3f\n",
            method, scale, seed, N, out$detection_mean, out$fidelity_mean))
