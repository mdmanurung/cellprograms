source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
source("/mnt/shared-workspace/cmp/bench/cand_v2.R")
cat("registered:", "s_subspace" %in% names(REG_V2), "s_subspace_ct95" %in% names(REG_V2), "\n")

run_one <- function(f, gate) {
  o <- readRDS(f)
  fit <- o$fit
  if (exists("get_profiles")) fit$profiles <- get_profiles(o)
  if (!is.null(o$profiles)) fit$profiles <- o$profiles
  t0 <- Sys.time()
  cl <- cand_subspace(fit, o$strata, 1L, gate = gate, n_perm = 199L)
  rt <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("%s | %s | gate=%s | %d rows | %d clusters | %.1fs\n",
              basename(f), o$truth$scenario, gate, nrow(cl), length(unique(cl$cluster)), rt))
  cl
}

f_mix <- "/mnt/shared-workspace/cmp/fits/K30/semi/cmp_mixed/s0001.rds"
f_nul <- "/mnt/shared-workspace/cmp/fits/K30/semi/null/s0001.rds"
f_par <- "/mnt/shared-workspace/cmp/fits/K30/param/cmp_mixed/s0001.rds"
for (f in c(f_mix, f_nul, f_par)) if (!file.exists(f)) cat("MISSING:", f, "\n")

cl1 <- run_one(f_mix, "maxT"); print(cl1)
cl2 <- run_one(f_mix, "ct95")
cl3 <- run_one(f_nul, "maxT")
cl4 <- run_one(f_nul, "ct95")
cl5 <- run_one(f_par, "maxT"); print(cl5)
