## arms.R -- one function per arm. Every arm returns a canonicalized cell_program_fit
## (plus a `use_strata` flag); run_arms() evaluates them all on one dataset (paired).
## Candidates are built only from existing arguments; nothing here changes R/.

INST <- "institute"
.args0 <- function() list(loading_prior = "point_laplace", var_type = 1, features = "variable",
                          n_variable_genes = getOption("ms.nvar", 2000L),
                          covariate_mode = "residualize", seed = 1L, max_factors = 10L)

.fit_joint <- function(x, covariates = INST, ...) {
  a <- modifyList(.args0(), list(...))
  canonicalize_programs(do.call(fit_celltype_programs, c(list(x, covariates = covariates), a)))
}

## fit each cell type separately (own K / features / S) and merge into one fit object
.fit_perct <- function(x, spec, covariates = INST) {
  fits <- lapply(setNames(x$cell_types, x$cell_types), function(ct) {
    xc <- as_cell_program_data(setNames(list(x$pseudobulk[[ct]]), ct),
                               sample_metadata = x$sample_metadata)
    a <- modifyList(.args0(), spec(ct))
    do.call(fit_celltype_programs, c(list(xc, covariates = covariates), a))
  })
  out <- x
  out$fits <- lapply(setNames(names(fits), names(fits)), function(ct) fits[[ct]]$fits[[ct]])
  out$preprocessing <- fits[[1]]$preprocessing
  out$preprocessing$feature_sets <- lapply(fits, function(f) f$preprocessing$feature_sets[[1]])
  out$provenance <- fits[[1]]$provenance
  class(out) <- c("cell_program_fit", "cell_program_data", "list")
  canonicalize_programs(out)
}

.resid_var <- function(fit, ct) {
  R <- fit$fits[[ct]]$Y_used - fit$scores[[ct]] %*% t(fit$loadings[[ct]])
  stats::median(apply(R, 2, stats::var))
}

arm_k_n <- function(x, ...) .fit_perct(x, function(ct) list(max_factors = as.integer(min(10L, floor(nrow(x$pseudobulk[[ct]]) / 10)))))
arm_k_20 <- function(x, ...) .fit_joint(x, max_factors = 20L)
arm_vt_12 <- function(x, ...) .fit_joint(x, var_type = c(1, 2))
arm_none <- function(x, ...) .fit_joint(x, covariates = NULL)
arm_genes_after <- function(x, ...) {
  nv <- getOption("ms.nvar", 2000L)
  .fit_perct(x, function(ct) {
    Y <- .residualize(x$pseudobulk[[ct]], x$sample_metadata, INST, ct)
    v <- apply(Y, 2, stats::var)
    list(features = setNames(list(names(sort(v, decreasing = TRUE))[seq_len(min(nv, length(v)))]), ct))
  })
}
## S_i = sqrt(c / n_i) with c set so median S^2 = 0.5 * median residual variance of the base fit
arm_w_S <- function(x, base_fit, ncell, ...) {
  .fit_perct(x, function(ct) {
    n_i <- ncell[[ct]][rownames(x$pseudobulk[[ct]])]
    S <- unname(sqrt(0.5 * .resid_var(base_fit, ct) * stats::median(n_i) / n_i))
    list(flash_control = list(S = S), var_type = 2)
  })
}
arm_k_stab <- function(x, base_fit, n_boot = getOption("ms.nboot", 10L), ...) {
  st <- assess_program_stability(base_fit, n_boot = n_boot, seed = 1, progress = FALSE)$stability$per_program
  keep <- st$program_id[!is.na(st$recovery_freq) & st$recovery_freq >= 0.7]
  fit <- base_fit
  for (ct in fit$cell_types) {
    cols <- intersect(colnames(fit$scores[[ct]]), keep)
    fit$scores[[ct]] <- fit$scores[[ct]][, cols, drop = FALSE]
    fit$loadings[[ct]] <- fit$loadings[[ct]][, cols, drop = FALSE]
  }
  fit
}

ARMS <- c("base", "k_n", "k_20", "k_stab", "genes_after", "vt_12", "w_S", "strata", "none", "none_strata")

.best_cor <- function(S, v) if (ncol(S) > 0L && length(v) > 2L) max(abs(suppressWarnings(stats::cor(S, v))), na.rm = TRUE) else NA_real_

## chance level of max |cor| over this arm's score columns: mean over 50 random z
.chance_cor <- function(S, n_draw = 50L) {
  if (ncol(S) == 0L) return(NA_real_)
  mean(replicate(n_draw, .best_cor(S, stats::rnorm(nrow(S)))))
}

.eval_fit <- function(fit, d, arm, use_strata, secs, n_perm) {
  tr <- d$truth
  ct_rows <- do.call(rbind, lapply(fit$cell_types, function(ct) {
    S <- fit$scores[[ct]]
    z <- tr$z[[ct]]
    data.frame(arm = arm, ct = ct, n = nrow(S), k = ncol(S),
               best_cor = if (is.null(z)) NA_real_ else .best_cor(S, z[rownames(S)]),
               chance = .chance_cor(S), fit_secs = secs, stringsAsFactors = FALSE)
  }))
  fail_msg <- ""
  ss <- tryCatch(sharing_spectrum(fit, space = "scores", n_perm = n_perm, seed = 1L,
                                  strata = if (use_strata) d$strata else NULL)$summary,
                 error = function(e) { fail_msg <<- conditionMessage(e); message("sharing failed (", arm, "): ", fail_msg); NULL })
  all_pairs <- utils::combn(fit$cell_types, 2L, function(p) paste(p, collapse = " vs "))
  pair_rows <- data.frame(arm = arm, pair = all_pairs, stringsAsFactors = FALSE)
  m <- if (is.null(ss)) NULL else match(pair_rows$pair, ss$pair)
  get <- function(col, default) if (is.null(m)) rep(default, nrow(pair_rows)) else ifelse(is.na(m), default, ss[[col]][m])
  ## tested = the sharing test produced a row for this pair. Untested is not a tested negative:
  ## `shared` stays FALSE (preregistered F unchanged), summarize.R gates on frac_untested.
  pair_rows$tested <- if (is.null(m)) rep(FALSE, nrow(pair_rows)) else !is.na(m)
  pair_rows$fail_reason <- if (!is.null(m)) ifelse(is.na(m), "pair absent (a cell type kept no program)", "") else fail_msg
  pair_rows$k_a <- get("k_a", NA_real_); pair_rows$k_b <- get("k_b", NA_real_)
  pair_rows$p_pair <- get("p_pair", NA_real_); pair_rows$q_pair <- get("q_pair", NA_real_)
  pair_rows$excess_frac <- get("excess_frac", NA_real_); pair_rows$z_pair <- get("z_pair", NA_real_)
  pair_rows$underpowered <- as.logical(get("underpowered", NA))   # NA = pair absent (a cell type kept 0 programs)
  pair_rows$shared <- as.logical(get("shared", FALSE)) & !is.na(pair_rows$underpowered)
  pair_rows$callable <- !is.na(pair_rows$underpowered) & !pair_rows$underpowered
  true_pair <- if (is.null(tr$pair) || tr$scenario == "S7") NA_character_ else paste(tr$pair, collapse = " vs ")
  pair_rows$is_true <- !is.na(true_pair) & pair_rows$pair == true_pair
  list(ct = ct_rows, pair = pair_rows)
}

run_arms <- function(d, arms = ARMS, n_perm = getOption("ms.nperm", 999L)) {
  x <- d$x
  out <- list(ct = list(), pair = list())
  cache <- new.env()
  timed <- function(expr) { t0 <- Sys.time(); v <- force(expr); list(v = v, s = as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
  base <- timed(.fit_joint(x)); cache$base <- base
  none <- NULL
  for (arm in arms) {
    r <- switch(arm,
      base = base,
      strata = base,
      k_n = timed(arm_k_n(x)),
      k_20 = timed(arm_k_20(x)),
      vt_12 = timed(arm_vt_12(x)),
      genes_after = timed(arm_genes_after(x)),
      w_S = timed(arm_w_S(x, base$v, d$ncell)),
      k_stab = { r0 <- timed(arm_k_stab(x, base$v)); r0$s <- r0$s + base$s; r0 },
      none = { if (is.null(none)) none <- timed(arm_none(x)); none },
      none_strata = { if (is.null(none)) none <- timed(arm_none(x)); none },
      stop("unknown arm ", arm))
    e <- .eval_fit(r$v, d, arm, arm %in% c("strata", "none_strata"), r$s, n_perm)
    out$ct[[arm]] <- e$ct; out$pair[[arm]] <- e$pair
  }
  list(ct = do.call(rbind, out$ct), pair = do.call(rbind, out$pair))
}
