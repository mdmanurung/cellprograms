source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
source("/mnt/shared-workspace/cmp/bench/cand_v2.R")

f <- "/mnt/shared-workspace/cmp/fits/K30/semi/cmp_mixed/s0002.rds"
cat("file size MB:", file.size(f) / 1e6, "\n")
t0 <- Sys.time(); o <- readRDS(f); cat("readRDS:", Sys.time() - t0, "\n")

fit <- o$fit
if (exists("get_profiles")) fit$profiles <- get_profiles(o)
t0 <- Sys.time(); cl <- cand_subspace(fit, o$strata, 1L, gate = "maxT", n_perm = 199L); cat("cand:", Sys.time() - t0, "\n")
t0 <- Sys.time(); ev <- eval_calls(cl, fit, o$truth, o$strata, "marginal"); cat("eval_calls:", Sys.time() - t0, "\n")
t0 <- Sys.time(); p <- get_profiles(o); cat("get_profiles:", Sys.time() - t0, "\n")

## time the full eval_method inner loop on 8 files, 8 cores
files <- file.path("/mnt/shared-workspace/cmp/fits/K30/semi/cmp_mixed", sprintf("s%04d.rds", 1:8))
t0 <- Sys.time()
res <- parallel::mclapply(files, function(f) {
  o <- readRDS(f); fit <- o$fit
  if (exists("get_profiles")) fit$profiles <- get_profiles(o)
  cl <- cand_subspace(fit, o$strata, 1L, gate = "maxT", n_perm = 199L)
  ev <- eval_calls(cl, fit, o$truth, o$strata, "marginal")
  nrow(ev$summary)
}, mc.cores = 8)
cat("8 files / 8 cores:", Sys.time() - t0, "\n")
