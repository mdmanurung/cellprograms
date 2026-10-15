log <- function(...) { cat(sprintf(...), "\n"); flush.console() }
source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
source("/mnt/shared-workspace/cmp/bench/cand_v2.R")
log("sourced")

f <- "/mnt/shared-workspace/cmp/fits/K30/semi/cmp_mixed/s0005.rds"
t0 <- Sys.time(); o <- readRDS(f); log("readRDS fit: %.1f", Sys.time() - t0)
t0 <- Sys.time(); pr <- get_profiles(o); log("get_profiles semi (cold): %.1f", Sys.time() - t0)
t0 <- Sys.time(); pr <- get_profiles(o); log("get_profiles semi (warm): %.2f", Sys.time() - t0)
fit <- o$fit; fit$profiles <- pr
t0 <- Sys.time(); cl <- cand_subspace(fit, o$strata, 1L, gate = "maxT", n_perm = 199L); log("cand: %.1f rows=%d", Sys.time() - t0, nrow(cl))
t0 <- Sys.time(); ev <- eval_calls(cl, fit, o$truth, o$strata, "marginal"); log("eval_calls: %.1f", Sys.time() - t0)
t0 <- Sys.time(); ev0 <- eval_calls(.empty_clusters(), fit, o$truth, o$strata, "marginal"); log("eval_calls empty: %.1f", Sys.time() - t0)

fp <- "/mnt/shared-workspace/cmp/fits/K30/param/cmp_mixed/s0002.rds"
t0 <- Sys.time(); op <- readRDS(fp); log("readRDS param fit: %.1f", Sys.time() - t0)
t0 <- Sys.time(); pp <- get_profiles(op); log("get_profiles param: %.2f", Sys.time() - t0)
fitp <- op$fit; fitp$profiles <- pp
t0 <- Sys.time(); clp <- cand_subspace(fitp, op$strata, 1L, gate = "maxT", n_perm = 199L); log("cand param: %.1f rows=%d", Sys.time() - t0, nrow(clp))
t0 <- Sys.time(); evp <- eval_calls(clp, fitp, op$truth, op$strata, "marginal"); log("eval_calls param: %.1f", Sys.time() - t0)
log("done")
