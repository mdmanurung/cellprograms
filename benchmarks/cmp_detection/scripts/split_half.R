## split_half.R -- split-half donor replication with programs re-learned in each half (Stephenson, institute+sex).
## Rscript split_half.R <n_splits> <cores>
source("/mnt/shared-workspace/cmp/bench/common.R")
for (f in c("coordination.R", "cmp-annotate.R")) sys.source(file.path(CMP_ROOT, "pkg/R", f), envir = globalenv())
suppressPackageStartupMessages(library(parallel))
a <- commandArgs(TRUE); ns <- as.integer(a[1]); cores <- as.integer(a[2])
b <- readRDS(file.path(CMP_ROOT, "data/stephenson_pb/semi_base.rds"))
full <- readRDS(file.path(CMP_ROOT, "real/stephenson_pb_sex_fit.rds"))
fres <- readRDS(file.path(CMP_ROOT, "real/stephenson_pb_sex_cmp.rds"))$flag
md <- full$sample_metadata
prof <- lapply(b$cts, function(z) log2(z$cpm_mean + 1))
jobs <- expand.grid(half = 1:2, split = seq_len(ns))
run <- function(i) {
  s <- jobs$split[i]; h <- jobs$half[i]
  set.seed(1000 + s)
  ids <- unlist(lapply(split(md$observation_id, md$institute), function(v) v[sample.int(length(v))][seq_len(ceiling(length(v) / 2))]))
  keep <- if (h == 1) ids else setdiff(md$observation_id, ids)
  Y <- lapply(b$cts, function(z) z$Y[intersect(rownames(z$Y), keep), , drop = FALSE])
  Y <- Y[vapply(Y, nrow, 1L) >= 20]
  m <- md[md$observation_id %in% keep, ]
  x <- as_cell_program_data(Y, sample_metadata = m)
  fit <- canonicalize_programs(fit_celltype_programs(x, loading_prior = "point_laplace", var_type = 1, features = "variable",
        n_variable_genes = 2000L, covariate_mode = "residualize", seed = 1L, max_factors = 10L, covariates = c("institute", "sex")))
  r <- find_coordinated_programs(fit, strata = "institute", ambient = "flag", profiles = prof, n_perm = 499)
  list(split = s, half = h, n = length(keep), fit_loadings = fit$loadings, cmps = r$cmps)
}
out <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(run(i), error = function(e) conditionMessage(e)), mc.cores = cores, mc.preschedule = FALSE)
tmp <- tempfile(fileext = ".rds"); saveRDS(out, tmp); system2("cp", c(tmp, file.path(CMP_ROOT, "real/split_half.rds")))
## match each full-data CMP to its best half-data CMP: mean over shared cell types of |cor| of member loadings
fc <- as.data.table(fres$cmps)
res <- rbindlist(lapply(out, function(o) {
  if (is.character(o)) return(NULL)
  hc <- as.data.table(o$cmps)
  rbindlist(lapply(unique(fc$cmp), function(i) {
    a <- fc[cmp == i]
    if (!nrow(hc)) return(data.table(split = o$split, half = o$half, cmp = i, best = NA, r = 0, n_ct = 0L))
    sc <- vapply(unique(hc$cmp), function(j) {
      bb <- hc[cmp == j]; cts <- intersect(a$cell_type, bb$cell_type); if (!length(cts)) return(0)
      rr <- vapply(cts, function(ct) {
        wa <- full$loadings[[ct]][, a[cell_type == ct]$program_id] * a[cell_type == ct]$sign
        wb <- o$fit_loadings[[ct]][, bb[cell_type == ct]$program_id] * bb[cell_type == ct]$sign
        g <- intersect(names(wa), names(wb)); if (length(g) < 50) return(NA_real_); stats::cor(wa[g], wb[g]) }, 1)
      abs(mean(rr, na.rm = TRUE)) * sum(!is.na(rr)) / nrow(a)   # penalize partial cell-type coverage
    }, 1)
    data.table(split = o$split, half = o$half, cmp = i, best = unique(hc$cmp)[which.max(sc)], r = max(sc), n_ct = nrow(a))
  }))
}))
suppressPackageStartupMessages(library(data.table))
agg <- res[, .(mean_score = round(mean(r), 3), frac_replicated = round(mean(r > 0.5), 2), n_halves = .N), by = cmp]
print(agg)
fwrite(res, "/tmp/sh.csv"); fwrite(agg, "/tmp/sha.csv")
system2("cp", c("/tmp/sh.csv", "/mnt/results/tusoskill/cmp_v1/real/stephenson_pb_sex/split_half_detail.csv"))
system2("cp", c("/tmp/sha.csv", "/mnt/results/tusoskill/cmp_v1/real/stephenson_pb_sex/split_half_summary.csv"))
