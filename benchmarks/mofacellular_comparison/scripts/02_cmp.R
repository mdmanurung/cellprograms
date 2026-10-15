## 02_cmp.R -- cross-cell-type coordination on a masked fit.
## Default settings (n_perm = 999, seed = 1, min_donors = 20), matching the
## published comparison. Run after 01_fit_masked.R (or 04_fit_excl.R).
##
## Usage: Rscript 02_cmp.R [fit_rds] [out_rds]
##   defaults: results/stephenson_K30_masked_fit.rds -> results/stephenson_K30_cmp.rds

suppressPackageStartupMessages(library(cellprograms))

a <- commandArgs(TRUE)
fit_rds <- if (length(a) >= 1) a[1] else file.path("results", "stephenson_K30_masked_fit.rds")
out_rds <- if (length(a) >= 2) a[2] else file.path("results", "stephenson_K30_cmp.rds")

f <- readRDS(fit_rds)
t0 <- Sys.time()
cmp <- find_coordinated_programs(f, n_perm = 999, seed = 1, min_donors = 20)
cat("CMP time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(cmp, out_rds)
cat("n CMPs:", length(unique(cmp$cmps$cmp)), "| edges:", nrow(cmp$edges), "\n")
