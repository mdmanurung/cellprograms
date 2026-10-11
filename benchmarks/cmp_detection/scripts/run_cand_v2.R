## run_cand_v2.R -- evaluate candidates (REG or REG_V2) on cached fits; v2 adds a context variable in marginal mode.
## Rscript run_cand_v2.R <cand1,...> <arm> <split: dev|test|heldout2|nullx> <cores> [bases=param,semi]
## env: CMP_SEEDS="1:15" restricts seeds; CMP_TAG suffix for output dir (default _v2); CMP_FITROOT (default fits)
source("/mnt/shared-workspace/cmp/bench/common.R")
source("/mnt/shared-workspace/cmp/bench/eval_cmp.R")
source("/mnt/shared-workspace/cmp/bench/config.R")
source("/mnt/shared-workspace/cmp/bench/registry.R")
source("/mnt/shared-workspace/cmp/bench/cand_v2.R")
ALL <- c(REG, REG_V2)

## context variable supplied to methods in marginal mode: true context (context scenario) or an independent dummy
## binary context (prevalence 0.4, standardized) drawn from the dataset seed -> tests false-call inflation of C2.
ctx_for <- function(o) {
  if (!is.null(o$context)) return(o$context)
  ids <- rownames(o$meta); if (is.null(ids)) ids <- o$fit$observation_ids
  set.seed(as.integer((as.numeric(o$truth$seed) * 7919 + 17) %% 2147483647))
  c0 <- as.numeric(stats::runif(length(ids)) < 0.4); if (stats::sd(c0) == 0) c0[1] <- 1 - c0[1]
  stats::setNames(as.numeric(scale(c0)), ids)
}

eval_method <- function(method_fun, files, cores = 16L, label = "m", seed = 1L) {
  res <- parallel::mclapply(files, function(f) {
    o <- readRDS(f)
    if (inherits(o$fit, "fit_error")) return(NULL)
    modes <- if (o$truth$scenario == "context") c("marginal", "within") else "marginal"
    lapply(modes, function(md) {
      fit <- if (md == "within") within_context(o$fit, o$context) else o$fit
      if (md == "within") { fit$scores_raw <- o$fit$scores; fit$context <- o$context } else fit$ctx_var <- ctx_for(o)
      if (exists("get_profiles")) fit$profiles <- get_profiles(o)
      if (!is.null(o$profiles)) fit$profiles <- o$profiles
      t0 <- Sys.time()
      st <- if (md == "within") within_strata(o$strata, o$context) else o$strata
      cl <- tryCatch(method_fun(fit, st, seed), error = function(e) structure(conditionMessage(e), class = "err"))
      rt <- as.numeric(Sys.time() - t0, units = "secs")
      if (!inherits(cl, "err") && !is.data.frame(cl)) cl <- cl$clusters
      if (inherits(cl, "err")) return(list(summary = data.frame(file = f, base = o$truth$base, scenario = o$truth$scenario,
                                                                 seed = o$truth$seed, mode = md, error = as.character(cl), runtime = rt)))
      ev <- eval_calls(cl, fit, o$truth, o$strata, md)
      list(summary = cbind(data.frame(file = f, base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed,
                                      mode = md, error = NA_character_, runtime = rt), ev$summary),
           signals = if (length(ev$f1)) data.frame(base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed,
                                                   mode = md, signal = names(ev$f1), f1 = ev$f1, oracle = ev$oracle) else NULL,
           calls = if (nrow(cl)) cbind(data.frame(base = o$truth$base, scenario = o$truth$scenario, seed = o$truth$seed, mode = md),
                                       cl, label = unname(ev$labels[as.character(cl$cluster)])) else NULL)
    })
  }, mc.cores = cores, mc.preschedule = FALSE)
  flat <- unlist(res, recursive = FALSE)
  bind <- function(k) data.table::rbindlist(lapply(flat, `[[`, k), fill = TRUE)
  list(summary = bind("summary"), signals = bind("signals"), calls = bind("calls"), method = label)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  cands <- strsplit(a[1], ",")[[1]]; arm <- a[2]; split <- a[3]; cores <- as.integer(a[4])
  bases <- if (length(a) > 4) strsplit(a[5], ",")[[1]] else c("param", "semi")
  seeds <- switch(split, dev = DEV_SEEDS, test = TEST_SEEDS, heldout2 = 1051:1100, nullx = 31:130)
  if (nzchar(Sys.getenv("CMP_SEEDS"))) seeds <- intersect(seeds, eval(parse(text = Sys.getenv("CMP_SEEDS"))))
  scens <- if (split == "nullx") "null" else CMP_SCEN
  fitroot <- file.path(CMP_ROOT, Sys.getenv("CMP_FITROOT", "fits"))
  files <- unlist(lapply(bases, function(b) unlist(lapply(scens, function(sc) file.path(fitroot, arm, b, sc, sprintf("s%04d.rds", seeds))))))
  files <- files[file.exists(files)]
  cat("files", length(files), "\n")
  outd <- file.path(CMP_ROOT, "evals", paste0(split, Sys.getenv("CMP_TAG", "_v2")), arm)
  dir.create(outd, recursive = TRUE, showWarnings = FALSE)
  for (cn in cands) {
    t0 <- Sys.time()
    ev <- eval_method(ALL[[cn]], files, cores = cores, label = cn)
    sc <- cbind(candidate = cn, arm = arm, split = split, n_files = length(files), scorecard(ev))
    tmp <- tempfile(fileext = ".rds"); saveRDS(ev, tmp); system2("cp", c(tmp, file.path(outd, paste0(cn, ".rds"))))
    tmpc <- tempfile(fileext = ".csv"); write.csv(sc, tmpc, row.names = FALSE); system2("cp", c(tmpc, file.path(outd, paste0(cn, "_score.csv"))))
    cat(cn, sprintf("%.0fs", as.numeric(Sys.time() - t0, units = "secs")), "\n")
    print(t(sc))
  }
}
