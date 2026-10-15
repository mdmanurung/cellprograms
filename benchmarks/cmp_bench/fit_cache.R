## fit_cache.R -- simulate + fit per-cell-type programs, cache fit essentials + truth.
## Rscript fit_cache.R <arm> <base> <scenarios,comma> <seed_from> <seed_to> <amp> <cores> [outroot]
## arms: main (repo M2 defaults, residualize on institute), noadj, K20, featall
source("/mnt/shared-workspace/cmp/bench/common.R")
suppressPackageStartupMessages(library(parallel))

ARMS <- list(
  main    = list(),
  noadj   = list(covariates = NULL),
  K20     = list(max_factors = 20L),
  featall = list(features = "all")
)
fit_args <- function(arm) modifyList(list(loading_prior = "point_laplace", var_type = 1, features = "variable",
                                          n_variable_genes = 2000L, covariate_mode = "residualize", seed = 1L,
                                          max_factors = 10L, covariates = "institute"), ARMS[[arm]], keep.null = TRUE)

## keep scores/loadings/ids + per-ct flash essentials (F_lfsr, L_psd, pve, elbo); drop data matrices
## (Y_used / pseudobulk are regenerable deterministically from the seed)
strip_fit <- function(fit) {
  out <- fit[c("cell_types", "observation_ids", "scores", "loadings", "canonicalization", "sample_metadata")]
  out$canonicalization <- lapply(fit$canonicalization, function(c) c[c("scale_multiplier", "sign_multiplier")])
  out$feature_sets <- fit$preprocessing$feature_sets
  out$settings <- fit$preprocessing[setdiff(names(fit$preprocessing), "feature_sets")]
  out$fits <- lapply(fit$fits, function(f) list(
    flash = f$flash[intersect(names(f$flash), c("n_factors", "pve", "elbo", "residuals_sd", "F_lfsr", "L_psd"))],
    genes = f$genes, fallback = f$fallback, elapsed_secs = f$elapsed_secs))
  class(out) <- class(fit)
  out
}

cache_one <- function(arm, base, sc, seed, amp, outroot) {
  f <- file.path(outroot, arm, base, sc, sprintf("s%04d.rds", seed))
  if (file.exists(f) && file.size(f) > 0) return(invisible("skip"))
  t0 <- Sys.time()
  d <- gen_cmp(base, sc, seed, amp)
  a <- fit_args(arm)
  fit <- tryCatch(canonicalize_programs(do.call(fit_celltype_programs, c(list(d$x), a))),
                  error = function(e) structure(conditionMessage(e), class = "fit_error"))
  rt <- as.numeric(Sys.time() - t0, units = "secs")
  obj <- list(fit = if (inherits(fit, "fit_error")) fit else strip_fit(fit), truth = d$truth, strata = d$strata,
              context = d$context, ncell = d$ncell, meta = d$x$sample_metadata, arm = arm, amp = amp, runtime = rt)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".rds"); saveRDS(obj, tmp)
  ok <- system2("cp", c(shQuote(tmp), shQuote(f))) == 0L; unlink(tmp)   # R file.copy gives 0-byte files on S3 FUSE
  if (!ok || file.size(f) == 0) stop("copy to cache failed: ", f)
  invisible(if (inherits(fit, "fit_error")) paste("ERR", fit) else sprintf("%.0fs", rt))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  arm <- a[1]; base <- a[2]; scs <- strsplit(a[3], ",")[[1]]; s0 <- as.integer(a[4]); s1 <- as.integer(a[5])
  amp <- as.numeric(a[6]); cores <- as.integer(a[7]); outroot <- if (length(a) > 7) a[8] else file.path(CMP_ROOT, "fits")
  jobs <- expand.grid(seed = s0:s1, sc = scs, stringsAsFactors = FALSE)
  res <- mclapply(seq_len(nrow(jobs)), function(i) {
    r <- tryCatch(cache_one(arm, base, jobs$sc[i], jobs$seed[i], amp, outroot), error = function(e) paste("ERR", conditionMessage(e)))
    cat(sprintf("%s %s %s s%d %s\n", arm, base, jobs$sc[i], jobs$seed[i], r)); r
  }, mc.cores = cores, mc.preschedule = FALSE)
  cat("finished", sum(grepl("ERR", unlist(res))), "errors\n")
}
