# Stage 3: mcRigor QC + geometry metrics for one membership file.
# Usage: Rscript benchmark_metrics.R <fixture_dir> <mem_dir> <scale> <seed> <method> <N> <out_csv>
# Pure R. mcRigor_DETECT consumes (single-cell Seurat obj, cell_membership df);
# ground-truth dubious at the metacell level: any member doublet or mixed types.

suppressPackageStartupMessages({
  library(Seurat)
  library(mcRigor)
  library(cluster)
})

args <- commandArgs(trailingOnly = TRUE)
fx_dir <- args[1]; mem_dir <- args[2]; scale <- args[3]; seed <- as.integer(args[4])
method <- args[5]; N <- as.integer(args[6]); out_csv <- args[7]
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

t_start <- Sys.time()

counts <- readRDS(file.path(fx_dir, paste0("counts_", scale, "_s", seed, ".rds")))
meta <- readRDS(file.path(fx_dir, paste0("meta_", scale, "_s", seed, ".rds")))
mem <- read.csv(file.path(mem_dir, sprintf("membership_%s_%s_s%d_N%d.csv",
                                           method, scale, seed, N)),
                stringsAsFactors = FALSE)
mem_vec <- mem$membership; names(mem_vec) <- mem$cell_id

obj <- CreateSeuratObject(counts = counts, assay = "RNA", min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
# mcRigor gets a RAW object: its buildmc aggregates counts itself and fails on
# a preprocessed object (empirically pinned: raw obj -> DETECT OK;
# ScaleData/PCA obj -> "RNA assay doesn't leave any cells" -> "Cannot find
# cells provided")
obj_raw <- CreateSeuratObject(counts = counts, assay = "RNA", min.cells = 0,
                              min.features = 0)

# --- mcRigor QC ---
# Calling convention pinned by direct test against the installed mcRigor HEAD:
# labels must be STRINGS that do NOT start with a digit. Pure integers fail
# DETECT's `!= ""` cell filter (coerces to NA); digit-leading strings trigger
# Seurat's `g`-prefix name mangling and crash buildmc ("No cell overlap").
# "MC<n>" labels pass both paths. Column named "1"; tgamma left at default
# (passing tgamma="1" re-triggers the integer-coercion crash in DETECT).
# NOTE: as.character() drops vector names in R 4.4 — build labels directly and
# reassign names explicitly (verbatim-verified construction).
lab <- paste0("MC", as.integer(factor(mem$membership,
                                      levels = sort(unique(mem$membership)))))
names(lab) <- mem$cell_id
cell_membership <- data.frame(X = lab, stringsAsFactors = FALSE,
                              check.names = FALSE, row.names = names(lab))
colnames(cell_membership) <- "1"
det <- mcRigor_DETECT(obj_singlecell = obj_raw, cell_membership = cell_membership,
                      Nrep = 50, gene_filter = 0.1, feature_use = 2000,
                      test_cutoff = 0.05, cor_method = "pearson",
                      aggregate_method = "sum", draw = FALSE, step_save = FALSE,
                      output_file = tempfile(fileext = ".rds"),
                      fields = character(0))
tab <- if (is.data.frame(det)) det else det$TabMC
est_dub <- as.character(tab$mcRigor)
if (is.null(est_dub) || all(is.na(est_dub))) {
  stop("mcRigor_DETECT returned no verdict column")
}
mc_names <- sub("^MC", "", rownames(tab))
code_to_lab <- sort(unique(lab))
code_idx <- suppressWarnings(as.integer(mc_names))
if (!any(is.na(code_idx)) && max(code_idx) <= length(code_to_lab)) {
  mc_names <- code_to_lab[code_idx]
}

# ground truth per metacell from planted labels: dubious if any member is an
# injected doublet or members span >1 planted type; purity = fraction of
# non-doublet cells of the dominant type.
# Compare per-cell LABELS (cell_lab), not raw membership values — mc_names are
# "MC<n>" strings while mem_vec holds raw codes/names for SuperCell/SEACells.
cell_lab <- lab
true_dub_vec <- vapply(mc_names, function(mc) {
  cells <- names(cell_lab)[cell_lab == mc]
  types <- meta$cell_type[match(cells, meta$cell_id)]
  dub <- meta$dubious[match(cells, meta$cell_id)]
  any(dub) || length(unique(types)) > 1
}, logical(1))
purity_true <- vapply(mc_names, function(mc) {
  cells <- names(cell_lab)[cell_lab == mc]
  types <- meta$cell_type[match(cells, meta$cell_id)]
  dub <- meta$dubious[match(cells, meta$cell_id)]
  tab_mc <- table(types[!dub])
  if (length(tab_mc) == 0) 0 else max(tab_mc) / length(cells)
}, numeric(1))
true_dub <- ifelse(true_dub_vec, "dubious", "trustworthy")

# aricode::ARI requires factors/integers — as.integer() on character labels
# produces NA inside sort_pairs
ari <- aricode::ARI(factor(true_dub), factor(est_dub))
n_true_dub <- sum(true_dub == "dubious")
recall <- if (n_true_dub > 0) sum(est_dub[true_dub == "dubious"] == "dubious") / n_true_dub else NA
n_est_dub <- sum(est_dub == "dubious")
precision <- if (n_est_dub > 0) sum(true_dub[est_dub == "dubious"] == "dubious") / n_est_dub else NA
ff <- if (!is.na(recall) && !is.na(precision) && (recall + precision) > 0)
  2 / (1 / recall + 1 / precision) else NA

# --- geometry (PCA space, pure R) ---
emb <- Embeddings(obj, "pca")[names(mem_vec), ]
d <- dist(emb)
sil <- silhouette(as.integer(factor(mem_vec)), d)
asw <- mean(sil[, 3])
# compactness: mean over PCs of within-metacell variance
comp <- vapply(sort(unique(mem_vec)), function(mc) {
  e <- emb[mem_vec == mc, , drop = FALSE]
  mean(apply(e, 2, var))
}, numeric(1))

# per-metacell table for the size-vs-purity confound figure
mc_sizes <- table(cell_lab)[mc_names]
pmc <- data.frame(
  method = method, scale = scale, seed = seed, target_N = N,
  metacell = mc_names,
  size = as.integer(mc_sizes),
  purity = purity_true,
  stringsAsFactors = FALSE
)
write.csv(pmc, sub("metrics_", "permc_", out_csv), row.names = FALSE)

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
out <- data.frame(
  method = method, scale = scale, seed = seed, target_N = N,
  n_metacells = length(unique(mem_vec)),
  purity_mean = mean(purity_true, na.rm = TRUE),
  purity_median = median(purity_true, na.rm = TRUE),
  ari_dub = ari, recall_dub = recall, precision_dub = precision, ff_dub = ff,
  asw = asw, compactness_mean = mean(comp), compactness_median = median(comp),
  runtime_min = elapsed,
  stringsAsFactors = FALSE
)
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("METRICS OK %s %s s%d N=%d: ARI=%.3f FF=%.3f ASW=%.3f\n",
            method, scale, seed, N, ari, ff, asw))
