## run_subspace.R -- multi-candidate eval driver reading fits from a LOCAL mirror (default
## /workspace/cmp_local/fits). Semantics identical to run_cand_v2.R (same modes, ctx variable, eval_calls,
## scorecard, output format and destination), but each fit file is read once and all requested candidates
## are evaluated on it -- avoids re-reading the slow S3-FUSE fit store per candidate.
## usage: Rscript run_subspace.R <cands,comma> <arm> <split: dev|test|heldout2|nullx> <cores> [bases=param,semi] [localroot]
## env: CMP_SEEDS, CMP_TAG (default _v2) as in run_cand_v2.R
source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
source("/mnt/shared-workspace/cmp/bench/cand_v2.R")
ALL <- c(REG, REG_V2)

ctx_for <- function(o) {
  if (!is.null(o$context)) return(o$context)
  ids <- rownames(o$meta); if (is.null(ids)) ids <- o$fit$observation_ids
  set.seed(as.integer((as.numeric(o$truth$seed) * 7919 + 17) %% 2147483647))
  c0 <- as.numeric(stats::runif(length(ids)) < 0.4); if (stats::sd(c0) == 0) c0[1] <- 1 - c0[1]
  stats::setNames(as.numeric(scale(c0)), ids)
}

a <- commandArgs(TRUE)
cands <- strsplit(a[1], ",")[[1]]; arm <- a[2]; split <- a[3]; cores <- as.integer(a[4])
bases <- if (length(a) > 4) strsplit(a[5], ",")[[1]] else c("param", "semi")
localroot <- if (length(a) > 5) a[6] else "/workspace/cmp_local/fits"
seeds <- switch(split, dev = DEV_SEEDS, test = TEST_SEEDS, heldout2 = 1051:1100, nullx = 31:130)
if (nzchar(Sys.getenv("CMP_SEEDS"))) seeds <- intersect(seeds, eval(parse(text = Sys.getenv("CMP_SEEDS"))))
scens <- if (split == "nullx") "null" else CMP_SCEN
files <- unlist(lapply(bases, function(b) unlist(lapply(scens, function(sc) file.path(localroot, arm, b, sc, sprintf("s%04d.rds", seeds))))))
files <- files[file.exists(files)]
cat("files", length(files), "\n"); flush.console()

## warm the semi profile cache in the parent so forked workers inherit it (avoids 16x slow FUSE reads)
if ("semi" %in% bases && exists("get_profiles")) invisible(get_profiles(list(truth = list(base = "semi"))))
cat("cache warm\n"); flush.console()

res <- parallel::mclapply(files, function(f) {
  o <- readRDS(f)
  if (inherits(o$fit, "fit_error")) return(NULL)
  modes <- if (o$truth$scenario == "context") c("marginal", "within") else "marginal"
  out <- list()
  for (md in modes) {
    fit <- if (md == "within") within_context(o$fit, o$context) else o$fit
    if (md == "within") { fit$scores_raw <- o$fit$scores; fit$context <- o$context } else fit$ctx_var <- ctx_for(o)
    if (exists("get_profiles")) fit$profiles <- get_profiles(o)
    if (!is.null(o$profiles)) fit$profiles <- o$profiles
    st <- if (md == "within") within_strata(o$strata, o$context) else o$strata
    for (cn in cands) {
      t0 <- Sys.time()
      cl <- tryCatch(ALL[[cn]](fit, st, 1L), error = function(e) structure(conditionMessage(e), class = "err"))
      rt <- as.numeric(Sys.time() - t0, units = "secs")
      if (!inherits(cl, "err") && !is.data.frame(cl)) cl <- cl$clusters
      if (inherits(cl, "err")) {
        out[[length(out) + 1]] <- list(cand = cn, summary = data.frame(file = f, base = o$truth$base, scenario = o$truth$scenario,
                                                                       seed = o$truth$seed, mode = md, error = as.character(cl), runtime = rt))
      } else {
        ev <- eval_calls(cl, fit, o$truth, o$strata, md)
        out[[length(out) + 1]] <- list(cand = cn,
          summary = cbind(data.frame(file = f, base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed,
                                     mode = md, error = NA_character_, runtime = rt), ev$summary),
          signals = if (length(ev$f1)) data.frame(base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed,
                                                  mode = md, signal = names(ev$f1), f1 = ev$f1, oracle = ev$oracle) else NULL,
          calls = if (nrow(cl)) cbind(data.frame(base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed, mode = md),
                                      cl, label = unname(ev$labels[as.character(cl$cluster)])) else NULL)
      }
    }
  }
  out
}, mc.cores = cores, mc.preschedule = FALSE)

flat <- unlist(res, recursive = FALSE)
outd <- file.path(CMP_ROOT, "evals", paste0(split, Sys.getenv("CMP_TAG", "_v2")), arm)
dir.create(outd, recursive = TRUE, showWarnings = FALSE)
for (cn in cands) {
  sub <- flat[vapply(flat, function(x) identical(x$cand, cn), TRUE)]
  bind <- function(k) data.table::rbindlist(lapply(sub, `[[`, k), fill = TRUE)
  ev <- list(summary = bind("summary"), signals = bind("signals"), calls = bind("calls"), method = cn)
  sc <- cbind(candidate = cn, arm = arm, split = split, n_files = length(files), scorecard(ev))
  tmp <- tempfile(fileext = ".rds"); saveRDS(ev, tmp); system2("cp", c(tmp, file.path(outd, paste0(cn, ".rds"))))
  tmpc <- tempfile(fileext = ".csv"); write.csv(sc, tmpc, row.names = FALSE); system2("cp", c(tmpc, file.path(outd, paste0(cn, "_score.csv"))))
  cat(cn, "done\n"); print(t(sc)); flush.console()
}
