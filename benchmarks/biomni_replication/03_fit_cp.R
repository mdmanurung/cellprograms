## 03_fit_cp.R <cohort> <K> <config: biomni|repo> <out_rds>
## L0 (no borrowing) per-cell-type EBMF with the fork's fit_celltype_programs, then canonicalize.
##  biomni: fork defaults (point_normal loading prior, normal score prior, var_type 2, backfit, nullcheck,
##          center=TRUE, scale=FALSE, all genes), max_factors = K (fork default K = 30).
##  repo  : point_laplace, var_type 1, backfit TRUE, top-2000 variable genes/CT (repo R/ebmf.R default), K (default 10).
a <- commandArgs(TRUE)
if (length(a) != 4) stop("usage: Rscript 03_fit_cp.R <cohort> <K> <biomni|repo> <out_rds>")
cohort <- a[1]; K <- as.integer(a[2]); config <- match.arg(a[3], c("biomni", "repo")); out <- a[4]
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
pkgload::load_all(file.path(here, "fork"), quiet = TRUE)

d <- file.path(Sys.getenv("BIOMNI_DATA", file.path(here, "data")), cohort)
cf <- Sys.glob(file.path(d, "mat", "ct_*.csv")); stopifnot(length(cf) > 0)
mats <- lapply(cf, function(f) as.matrix(read.csv(f, row.names = 1, check.names = FALSE)))
names(mats) <- sub("\\.csv$", "", sub("^ct_", "", basename(cf)))
meta <- read.csv(file.path(d, "donor_meta.csv"), row.names = 1, check.names = FALSE)
meta <- data.frame(observation_id = rownames(meta), meta, check.names = FALSE)

x <- as_cell_program_data(mats, sample_metadata = meta,
                          features = if (config == "repo") "variable" else "all")
args <- if (config == "biomni") {
  list(loading_prior = "point_normal", var_type = 2L, backfit = TRUE, nullcheck = TRUE)
} else {
  list(loading_prior = "point_laplace", var_type = 1L, backfit = TRUE, nullcheck = TRUE, n_variable = 2000L)
}
t0 <- Sys.time()
fit <- do.call(fit_celltype_programs, c(list(x, max_factors = K, seed = 1L), args))
fit <- canonicalize_programs(fit)
fit$provenance$config <- list(name = config, K = K, cohort = cohort, args = args,
                              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
print(fit)
saveRDS(fit, out)
cat(sprintf("saved %s (%.1f min)\n", out, fit$provenance$config$minutes))
