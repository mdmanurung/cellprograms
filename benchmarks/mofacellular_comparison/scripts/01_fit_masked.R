## 01_fit_masked.R -- artifact-masked K30 refit of the Stephenson pseudobulk.
## Produces the default-setting fit used by the comparison vignette.
##
## Input: a pseudobulk object with per-cell-type donor x gene count matrices and
## donor metadata (see README.md for how to build it from ArrayExpress
## E-MTAB-8731). Set PB_DIR to the directory holding stephenson_pb/semi_base.rds.
##
## Output: <OUT_DIR>/stephenson_K30_masked_fit.rds
##
## Recipe: point-Laplace loading prior, K = 30, var_type = 1, backfit +
## nullcheck, covariates institute + sex residualized, 2000 HVGs, seed 1,
## exclude_genes = artifact_genes(universe, sets = c("vdj", "genotype_marker"))
## (the 274-gene default blacklist). Exclusion happens before HVG selection.

suppressPackageStartupMessages(library(cellprograms))

PB_DIR  <- Sys.getenv("PB_DIR",  "data")          # holds stephenson_pb/semi_base.rds
OUT_DIR <- Sys.getenv("OUT_DIR", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

b <- readRDS(file.path(PB_DIR, "stephenson_pb", "semi_base.rds"))
Y <- lapply(b$cts, `[[`, "Y")
ids <- sort(unique(unlist(lapply(Y, rownames))))
dm <- b$donor_meta; dm <- dm[match(ids, dm$donor_id), ]
grp <- dm$Status_on_day_collection_summary
sev <- c(Healthy = 0, Asymptomatic = 1, Mild = 2, Moderate = 3,
         Severe = 4, Critical = 5)[grp]
grp2 <- ifelse(grepl("^LPS", grp), "LPS",
               ifelse(grp == "Non_covid", "Non_covid",
                      ifelse(grp == "Healthy", "Healthy", "COVID")))
meta <- data.frame(observation_id = ids, institute = factor(dm$site),
                   group = grp2, group_fine = grp, severity = unname(sev),
                   sex = dm$sex, stringsAsFactors = FALSE)

all_genes <- unique(unlist(lapply(Y, colnames)))
excl <- artifact_genes(all_genes, sets = c("vdj", "genotype_marker"))
cat("gene universe:", length(all_genes), "| excluded:", length(excl), "\n")

x <- as_cell_program_data(Y, sample_metadata = meta)
t0 <- Sys.time()
fit <- canonicalize_programs(fit_celltype_programs(
  x, loading_prior = "point_laplace", var_type = 1, features = "variable",
  n_variable_genes = 2000L, exclude_genes = excl, max_factors = 30L,
  covariates = c("institute", "sex"), covariate_mode = "residualize", seed = 1L))
cat("fit time:", format(Sys.time() - t0), "\n")
fit$pseudobulk <- Y
saveRDS(fit, file.path(OUT_DIR, "stephenson_K30_masked_fit.rds"))
cat("programs:", sum(vapply(fit$loadings, ncol, 1L)), "\n")
