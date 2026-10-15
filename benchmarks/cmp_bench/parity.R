## parity: package find_coordinated_programs() vs benchmark n_agg_noanc on dev fits
source("/mnt/shared-workspace/cmp/bench/common.R"); source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R"); source("/mnt/shared-workspace/cmp/bench/registry.R")
for (f in list.files("/mnt/shared-workspace/cmp/pkg/R", full.names = TRUE)) source(f)
fs <- c(sprintf("/mnt/shared-workspace/cmp/fits/main/%s/%s/s%04d.rds", rep(c("param","semi"), each = 4),
        rep(c("cmp_distinct","weak","ambient","null"), 2), c(1,2,3,4,5,6,7,8)),
        "/mnt/shared-workspace/cmp/fits/main/param/context/s0002.rds")
for (fn in fs) {
  o <- readRDS(fn); f <- o$fit; f$profiles <- get_profiles(o)
  f$observation_ids <- sort(unique(unlist(lapply(f$scores, rownames))))
  a <- REG$n_agg_noanc(f, o$strata, 1)
  b <- find_coordinated_programs(f, strata = o$strata, profiles = lapply(f$profiles, function(p) log2(p + 1)), seed = 1)
  ka <- sort(paste(a$cell_type, a$program_id)); kb <- sort(paste(b$cmps$cell_type, b$cmps$program_id))
  pa <- sort(vapply(split(paste(a$cell_type, a$program_id), a$cluster), function(x) paste(sort(x), collapse = ";"), ""))
  pb <- if (nrow(b$cmps)) sort(vapply(split(paste(b$cmps$cell_type, b$cmps$program_id), b$cmps$cmp), function(x) paste(sort(x), collapse = ";"), "")) else character()
  cat(basename(dirname(fn)), o$truth$base, "bench", nrow(a), "pkg", nrow(b$cmps), "same_nodes", identical(ka, kb), "same_partition", identical(unname(pa), unname(pb)), "\n")
  if (o$truth$scenario == "context") {
    w <- find_coordinated_programs(f, strata = o$strata, context = o$context, mode = "within", profiles = lapply(f$profiles, function(p) log2(p + 1)))
    cat("  within mode CMPs:", length(unique(w$cmps$cmp)), "\n"); print(w)
  }
}
