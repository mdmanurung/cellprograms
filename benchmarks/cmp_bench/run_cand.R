## run_cand.R -- evaluate registry candidates on cached fits.
## Rscript run_cand.R <cand1,cand2,...> <arm> <split: dev|test> <cores> [bases=param,semi]
source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
a <- commandArgs(TRUE)
cands <- strsplit(a[1], ",")[[1]]; arm <- a[2]; split <- a[3]; cores <- as.integer(a[4])
bases <- if (length(a) > 4) strsplit(a[5], ",")[[1]] else c("param", "semi")
seeds <- switch(split, dev = DEV_SEEDS, test = TEST_SEEDS, nullx = 31:130)
scens <- if (split == "nullx") "null" else CMP_SCEN
files <- unlist(lapply(bases, function(b) unlist(lapply(scens, function(sc)
  file.path(CMP_ROOT, "fits", arm, b, sc, sprintf("s%04d.rds", seeds))))))
files <- files[file.exists(files)]
cat("files", length(files), "\n")
outd <- file.path(CMP_ROOT, "evals", paste0(split, Sys.getenv("CMP_TAG", "")), arm)   # CMP_TAG: new dir per benchmark version (S3 overwrite across machines fails)
dir.create(outd, recursive = TRUE, showWarnings = FALSE)
for (cn in cands) {
  t0 <- Sys.time()
  ev <- eval_method(REG[[cn]], files, cores = cores, label = cn)
  sc <- cbind(candidate = cn, arm = arm, split = split, n_files = length(files), scorecard(ev))
  tmp <- tempfile(fileext = ".rds"); saveRDS(ev, tmp); system2("cp", c(tmp, file.path(outd, paste0(cn, ".rds"))))
  tmpc <- tempfile(fileext = ".csv"); write.csv(sc, tmpc, row.names = FALSE); system2("cp", c(tmpc, file.path(outd, paste0(cn, "_score.csv"))))
  cat(cn, sprintf("%.0fs", as.numeric(Sys.time() - t0, units = "secs")), "\n")
  print(t(sc))
}
