## 07_adjusted_sharing.R <cohort> <covariates: none|col1,col2> <out_prefix>
## Fit the REPO fit_celltype_programs (K=10 defaults, top-2000 variable genes) with optional covariates=,
## then calibrated sharing spectrum (scores space, 1999 perms, BH over callable pairs).
## Fits both arms with the same function so adjusted vs unadjusted differ only in covariates=.
a <- commandArgs(TRUE); stopifnot(length(a) == 3)
cohort <- a[1]; cov <- if (a[2] == "none") NULL else strsplit(a[2], ",")[[1]]; out <- a[3]
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
library(stats); library(utils)
for (f in list.files(file.path(here, "..", "..", "R"), full.names = TRUE)) source(f)

d <- file.path(here, "data", cohort)
cf <- Sys.glob(file.path(d, "mat", "ct_*.csv"))
mats <- lapply(cf, function(f) as.matrix(read.csv(f, row.names = 1, check.names = FALSE)))
names(mats) <- sub("\\.csv$", "", sub("^ct_", "", basename(cf)))
meta <- read.csv(file.path(d, "donor_meta.csv"), row.names = 1, check.names = FALSE)
meta <- data.frame(observation_id = rownames(meta), meta, check.names = FALSE)
## blank strings are not NA to read.csv; give them an explicit level (e.g. Outcome is blank for the Sepsis donors)
for (cv in intersect(cov, names(meta))) if (is.character(meta[[cv]])) meta[[cv]][meta[[cv]] == ""] <- "not_recorded"
x <- as_cell_program_data(mats, sample_metadata = meta)

t0 <- Sys.time()
fit <- canonicalize_programs(fit_celltype_programs(x, features = "variable", covariates = cov, seed = 1))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(fit, paste0(out, "_fit.rds"))
ss <- suppressWarnings(sharing_spectrum(fit, space = "scores", n_perm = 1999, seed = 1))$summary
write.csv(ss, paste0(out, "_pairs.csv"), row.names = FALSE)
cat(sprintf("covariates=%s fit %.1f min; pairs %d underpowered %d shared %d\n",
            a[2], mins, nrow(ss), sum(ss$underpowered), sum(ss$shared)))
