source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
source("/mnt/shared-workspace/cmp/bench/cand_v2.R")
f <- "/mnt/shared-workspace/cmp/fits/K30/semi/cmp_mixed/s0002.rds"
t0 <- Sys.time(); o <- readRDS(f); cat("readRDS:", round(Sys.time() - t0, 2), "\n"); flush.console()
fit <- o$fit
t0 <- Sys.time(); fit$profiles <- get_profiles(o); cat("get_profiles:", round(Sys.time() - t0, 2), "\n"); flush.console()
t0 <- Sys.time(); cl <- cand_subspace(fit, o$strata, 1L, gate = "maxT", n_perm = 199L); cat("cand:", round(Sys.time() - t0, 2), "rows:", nrow(cl), "\n"); flush.console()
t0 <- Sys.time(); ev <- eval_calls(cl, fit, o$truth, o$strata, "marginal"); cat("eval_calls:", round(Sys.time() - t0, 2), "\n"); flush.console()
t0 <- Sys.time(); ev0 <- eval_calls(.empty_clusters(), fit, o$truth, o$strata, "marginal"); cat("eval_calls(empty):", round(Sys.time() - t0, 2), "\n"); flush.console()
