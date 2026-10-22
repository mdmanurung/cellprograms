## test_eb_cp.R — eb_cp gate (P0-2 milestone): sanity recovery + null calibration.
## Also smoke-tests the full stage-1 -> block_E path at tiny scale (flashier).
## run: R_LIBS_USER=/workspace/.Rlib Rscript benchmarks/coordination/test_eb_cp.R
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); dirname(normalizePath(f)) })
root <- normalizePath(file.path(here, "..", ".."))
library(stats)
for (f in c("data-model.R", "utils.R", "ebmf.R", "principal-angles.R", "block-joint.R", "eb-cp.R")) source(file.path(root, "R", f))
for (f in c("sim_coord.R", "methods.R", "eb_cp.R", "scoring.R")) source(file.path(here, f))
chk <- function(ok, msg) if (!isTRUE(ok)) stop("FAILED: ", msg) else cat("ok:", msg, "\n")
options(ms.G = 1200L)

## 1. sanity: eb_cp recovers an act+genes R100 component
d <- gen_coord("act+genes", amp = 0.75, seed = 21, breadth = "R100")
mats <- adjusted_matrices(d$x, "institute")
t0 <- Sys.time()
cp <- eb_cp_fit(mats, Kmax = 12L, seed = 1)
cat(sprintf("eb_cp: %d components, %.1fs, frac_unexplained %.3f\n",
            cp$n_components, as.numeric(difftime(Sys.time(), t0, units = "secs")), cp$frac_unexplained))
cc <- eb_cp_clusters(cp)
s <- score_method(list(method = "eb_cp", clusters = cc$clusters, activity = cc$activity), d$truth)
chk(s$program_level$found[1], "eb_cp finds R100 (act+genes, amp 0.75)")
chk(s$program_level$recall[1] >= 0.8, "eb_cp R100 membership recall >= 0.8")

## 2. null calibration: 20 null datasets -> coordination calls in at most 1
n_null <- 20L; bad <- 0L
for (i in seq_len(n_null)) {
  d <- gen_coord("null", amp = 0.5, seed = 100 + i)
  cp <- eb_cp_fit(adjusted_matrices(d$x, "institute"), Kmax = 12L, seed = i)
  cc <- eb_cp_clusters(cp)
  if (nrow(cc$clusters) > 0L) bad <- bad + 1L
}
cat(sprintf("null calibration: %d/%d datasets with >=1 call\n", bad, n_null))
chk(bad <= 1L, "eb_cp null: calls in <= 1/20 datasets")

## 3. stage-1 + block_E smoke at tiny scale
d <- gen_coord("act+genes", amp = 0.75, seed = 31, breadth = "R100")
t0 <- Sys.time()
fit <- fit_celltype_programs(d$x, max_factors = 8L, covariates = "institute",
                             covariate_mode = "residualize", seed = 31)
fit <- canonicalize_programs(fit)
cat(sprintf("stage-1: %.1fs, programs/ct: %s\n", as.numeric(difftime(Sys.time(), t0, units = "secs")),
            paste(vapply(fit$scores, ncol, 1L), collapse = ",")))
t0 <- Sys.time()
be <- method_block_e(fit, max_factors = 12L, seed = 1)
cat(sprintf("block_E: %.1fs, %d calls\n", as.numeric(difftime(Sys.time(), t0, units = "secs")),
            length(unique(be$clusters$cluster))))
s <- score_method(be, d$truth)
chk(s$program_level$found[1], "block_E finds R100 (act+genes, amp 0.75)")

## 4. eb_cp via method wrapper on the same fit
t0 <- Sys.time()
ec <- method_eb_cp(fit, strata = d$strata, Kmax = 12L, seed = 1)
cat(sprintf("eb_cp (wrapper): %.1fs, %d calls\n", as.numeric(difftime(Sys.time(), t0, units = "secs")),
            length(unique(ec$clusters$cluster))))
s <- score_method(ec, d$truth)
chk(s$program_level$found[1], "eb_cp (wrapper on Y_used) finds R100")
cat("eb_cp gate checks passed\n")
