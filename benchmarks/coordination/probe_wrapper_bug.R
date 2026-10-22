## probe_wrapper_bug.R — isolate why method_eb_cp (Y_used path) finds 0 calls
## on the seed-31 dataset while eb_cp_fit(adjusted_matrices(...)) finds R100.
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/probe_wrapper_bug.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R")) source(file.path(here, f))
options(ms.G = 1200L)

d <- gen_coord("act+genes", amp = 0.75, seed = 31, breadth = "R100")
mats_adj <- adjusted_matrices(d$x, "institute")

fit <- fit_celltype_programs(d$x, max_factors = 8L, covariates = "institute",
                             covariate_mode = "residualize", seed = 31)
fit <- canonicalize_programs(fit)
mats_y <- lapply(fit$fits[fit$cell_types], function(f) f$Y_used)

cat("== input comparison (Y_used vs adjusted_matrices) ==\n")
for (ct in names(mats_adj)) {
  a <- mats_adj[[ct]]; b <- mats_y[[ct]]
  cg <- intersect(colnames(a), colnames(b)); cr <- intersect(rownames(a), rownames(b))
  mad <- if (length(cg) && length(cr)) max(abs(a[cr, cg] - b[cr, cg])) else NA
  cat(sprintf("%s: adj %dx%d | yused %dx%d | common donors %d genes %d | max|diff| %.3g\n",
              ct, nrow(a), ncol(a), nrow(b), ncol(b), length(cr), length(cg), mad))
}

run_case <- function(label, mats, strata) {
  t0 <- Sys.time()
  cp <- eb_cp_fit(mats, Kmax = 12L, seed = 1L, strata = strata)
  cc <- eb_cp_clusters(cp)
  cat(sprintf("%s: n_comp=%d, calls=%d, %.1fs\n   trace_z: %s\n",
              label, cp$n_components, length(unique(cc$clusters$cluster)),
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              paste(sprintf("%.2f", cp$trace_z), collapse = " ")))
  if (cp$n_components > 0L) {
    w <- cp$W[, 1]
    cat(sprintf("   CP01 |w|: %s\n", paste(sprintf("%.2f", w), collapse = " ")))
  }
  invisible(cp)
}

run_case("Y_used   + strata", mats_y, d$strata)
run_case("Y_used   + NULL  ", mats_y, NULL)
run_case("adjusted + strata", mats_adj, d$strata)
run_case("adjusted + NULL  ", mats_adj, NULL)
