## 05_sim_pa_vs_pmd.R -- principal angles (PA) vs multi-view PMD on simulated data.
## (v2: see "v2 changes" below)
##
## Usage: Rscript 05_sim_pa_vs_pmd.R <seed> <outfile.csv> [--smoke]
##   --smoke : scenarios S2, S4, S8, S9 only.
##
## Scenarios (5 cell types -> 10 pairs): S1-S7 from fork/R/simulate.R; S8 "confounded"
## = S4 + a strong binary donor-level nuisance; S9 = S3 + the same nuisance (added
## because S4's truth is all-pairs-shared, so S8 alone cannot show nuisance-induced
## false positives). Truth: S1,S5,S7 = 0 shared pairs (10 unshared); S2,S4,S6,S8 = 10
## shared; S3,S9 = 1 shared pair (CD8-NK) and 9 unshared. Only S3/S9 have both classes,
## so only they get AUROC / TPR@FPR. The nuisance is ONE donor vector (identical in
## every cell type, indexed by observation ID) acting through cell-type-specific
## sparse gene loadings, added to every cell type's matrix before stage-1 fitting.
##
## Arms: PA_scores, PA_loadings, PMD_raw; on S8/S9 only, nuisance-fix arms:
##   PMD_resid, PA_scores_resid   -- residualize scores on the ORACLE nuisance vector;
##   PMD_estnuis, PA_scores_estnuis -- residualize on an ESTIMATED nuisance (no truth).
##
## ESTIMATED-NUISANCE RULE (est_nuisance()). Using only the stage-1 scores:
##   1. for each cell type, PC1 of its standardized score matrix (prcomp, scale=TRUE;
##      a 1-factor cell type contributes that factor);
##   2. sign-align each cell type's PC1 to the first cell type's (sign of the
##      correlation over shared donors);
##   3. per donor, average the aligned PC1 scores over the cell types that have the
##      donor -> the "consensus donor axis" (one vector);
##   4. regress that axis out of every cell type's score columns (lm residuals).
## Caveat: the axis is whatever dominates each cell type's scores; if a true shared
## program dominates (e.g. S4's global program), it is removed too. The metric
## `nuis_est_cor` records |cor(estimated axis, true nuisance)| for the estnuis arms.
##
## TRUTH. Pair (a,b) is truly shared iff some program in $truth$structure with
## type != "private" lists BOTH a and b in its `cell_types` ("A+B+..." string).
## The true shared activity of the pair is $truth$activities[[prog]]$z for those
## programs (indexed by obs ID "obsN").
##
## PAIR CALLS AND SCORES.
##  PA : shared iff principal_angles()$n_shared_05 >= 1 (p < 0.05 AND cosine above
##       the random-subspace reference). score_cont = max cosine;
##       score_sig = score_cont if called else 0.
##  PMD: shared iff some MCP has perm p < 0.05 AND pair_cors[a,b] > 0.5 in it.
##       score_cont = max pair_cors[a,b] over ALL MCPs (-1 if the pair is in none);
##       score_sig  = max pair_cors[a,b] over significant MCPs (0 if none), i.e. the
##       v1 score, which ties every non-significant pair at 0 and understates AUROC.
##
## v2 changes: (1) a failed PA run is NA (not "not shared") and counted in `failed`;
## (2) AUROC reported with both scorings (auroc_sig, auroc_cont; for PA the v1
## `auroc` equals auroc_cont); (3) TPR at FPR <= 0.05 (tpr05_sig, tpr05_cont);
## (4) estimated-nuisance arms.
##
## METRICS (long CSV: seed, scenario, missing_frac, arm, metric, value; NA dropped):
##  auroc_sig, auroc_cont, tpr05_sig, tpr05_cont (over pairs; S3/S9 only), tpr, fpr
##  (pairs, over non-failed pairs), failed (fraction of pairs whose run failed; 1 if
##  the whole arm errored), n_false_shared (S1,S5,S7: number of called pairs, all of
##  which are false), state_cor (mean over truly shared pairs of the best |cor|
##  between a pair-level state and the true shared activity; PA: mean of paired
##  principal vectors i, MCP: mean of the pair's two views' MCP scores), seconds.
##  arm "stage1" holds fit seconds.

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args
args <- setdiff(args, "--smoke")
seed <- as.integer(args[1]); outfile <- args[2]
if (is.na(seed) || is.na(outfile)) stop("usage: Rscript 05_sim_pa_vs_pmd.R <seed> <outfile> [--smoke]")

here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
suppressPackageStartupMessages(pkgload::load_all(file.path(here, "fork"), quiet = TRUE))

N <- 100L; G <- 500L; N_PERM_PA <- 200L; N_PERM_PMD <- 50L
scen <- list(S1 = c(1, FALSE), S2 = c(2, FALSE), S3 = c(3, FALSE), S4 = c(4, FALSE),
             S5 = c(5, FALSE), S6 = c(6, FALSE), S7 = c(7, FALSE),
             S8 = c(4, TRUE), S9 = c(3, TRUE))
all_private <- c("S1", "S5", "S7")
if (smoke) scen <- scen[c("S2", "S4", "S8", "S9")]
if (nzchar(Sys.getenv("SIM_SCEN"))) scen <- scen[strsplit(Sys.getenv("SIM_SCEN"), ",")[[1]]]  # debug override

## ---- truth -----------------------------------------------------------------
truth_pairs <- function(truth, cts) {
  pr <- t(utils::combn(cts, 2L)); colnames(pr) <- c("a", "b")
  pr <- as.data.frame(pr, stringsAsFactors = FALSE)
  st <- truth$structure[truth$structure$type != "private", , drop = FALSE]
  pr$progs <- lapply(seq_len(nrow(pr)), function(i) {
    st$program[vapply(strsplit(st$cell_types, "+", fixed = TRUE),
                      function(v) all(c(pr$a[i], pr$b[i]) %in% v), logical(1L))]
  })
  pr$truth <- lengths(pr$progs) > 0L
  pr
}

## ---- nuisance --------------------------------------------------------------
add_nuisance <- function(x, seed, load_sd = 6, sparsity = 0.05) {
  set.seed(seed)
  ids <- sort(unique(unlist(lapply(x$matrices, rownames))))
  nuis <- stats::setNames(as.numeric(sample(rep(0:1, length.out = length(ids)))), ids)
  for (ct in names(x$matrices)) {
    Y <- x$matrices[[ct]]
    w <- rep(0, ncol(Y)); idx <- sample.int(ncol(Y), round(ncol(Y) * sparsity))
    w[idx] <- stats::rnorm(length(idx), 0, load_sd)
    x$matrices[[ct]] <- Y + nuis[rownames(Y)] %o% w
  }
  list(x = x, nuis = nuis)
}

residualize <- function(fit, nuis) {
  for (ct in names(fit$scores)) {
    Z <- fit$scores[[ct]]
    if (ncol(Z) == 0L) next
    fit$scores[[ct]] <- stats::resid(stats::lm(Z ~ nuis[rownames(Z)]))
  }
  fit
}

est_nuisance <- function(fit) {            # rule documented in the header
  pcs <- lapply(fit$scores, function(Z) {
    if (ncol(Z) == 0L) return(NULL)
    Z <- Z[, apply(Z, 2L, stats::sd) > 0, drop = FALSE]
    if (ncol(Z) == 0L) return(NULL)
    if (ncol(Z) == 1L) return(stats::setNames(as.numeric(scale(Z)), rownames(Z)))
    stats::setNames(stats::prcomp(Z, scale. = TRUE)$x[, 1L], rownames(Z))
  })
  pcs <- Filter(Negate(is.null), pcs)
  ref <- pcs[[1L]]
  pcs <- lapply(pcs, function(v) {
    k <- intersect(names(v), names(ref))
    r <- if (length(k) > 2L) stats::cor(v[k], ref[k]) else NA
    if (isTRUE(r < 0)) -v else v
  })
  ids <- sort(unique(unlist(lapply(pcs, names))))
  M <- vapply(pcs, function(v) v[ids], numeric(length(ids)))
  stats::setNames(rowMeans(M, na.rm = TRUE), ids)
}

## ---- arm runners: each returns data.frame(a, b, called, score_sig, score_cont, state_cor) ----
cor_best <- function(states, z) {          # states: matrix (obs x k) or vector
  states <- as.matrix(states)
  max(abs(stats::cor(states, z[rownames(states)], use = "pairwise.complete.obs")), na.rm = TRUE)
}

run_pa <- function(fit, pairs, space, tr) {
  res <- lapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs$a[i]; b <- pairs$b[i]
    pa <- tryCatch(suppressWarnings(principal_angles(fit, a, b, space = space,
                                                     n_perm = N_PERM_PA, seed = seed)),
                   error = function(e) NULL)
    if (is.null(pa)) return(c(called = NA, score_sig = NA, score_cont = NA, state_cor = NA))  # failure = NA, not "not shared"
    sc <- NA_real_
    if (pairs$truth[i]) {
      ## average the paired principal vectors across the pair (sign-aligned);
      ## in loadings space the vectors live in gene space -> no sample-level state
      if (space == "scores") {
        pa_a <- pa$principal_vectors$a; pa_b <- pa$principal_vectors$b
        sgn <- sign(colSums(pa_a * pa_b)); sgn[sgn == 0] <- 1
        st <- (pa_a + sweep(pa_b, 2L, sgn, `*`)) / 2
        rownames(st) <- intersect(rownames(fit$scores[[a]]), rownames(fit$scores[[b]]))
        sc <- max(vapply(pairs$progs[[i]], function(p) cor_best(st, tr$z[[p]]), numeric(1L)))
      }
    }
    called <- as.numeric(pa$n_shared_05 >= 1)
    c(called = called, score_sig = if (called == 1) pa$cosines[1] else 0,
      score_cont = pa$cosines[1], state_cor = sc)
  })
  cbind(pairs[, c("a", "b")], as.data.frame(do.call(rbind, res)))
}

run_pmd <- function(fit, pairs, tr) {
  fit <- suppressWarnings(integrate_programs(fit, method = "pmd", n_states = 3L,
                                             n_perm = N_PERM_PMD, seed = seed, verbose = 0L))
  mcps <- fit$integration$mcps
  res <- lapply(seq_len(nrow(pairs)), function(i) {
    a <- pairs$a[i]; b <- pairs$b[i]
    r <- vapply(mcps, function(m) {
      pc <- if (all(c(a, b) %in% rownames(m$pair_cors))) m$pair_cors[a, b] else NA_real_
      c(sig = as.numeric(isTRUE(m$pvalue < 0.05) && isTRUE(pc > 0.5)),
        sc = if (isTRUE(m$pvalue < 0.05) && !is.na(pc)) pc else NA_real_,
        pc = pc)
    }, numeric(3L))
    r <- matrix(r, nrow = 3L)
    sc <- NA_real_
    if (pairs$truth[i] && length(mcps) > 0L) {
      st <- vapply(mcps, function(m) rowMeans(m$state_by_celltype[, c(a, b), drop = FALSE],
                                              na.rm = TRUE), numeric(length(mcps[[1]]$state)))
      sc <- max(vapply(pairs$progs[[i]], function(p) cor_best(st, tr$z[[p]]), numeric(1L)))
    }
    c(called = as.numeric(any(r[1, ] == 1)),
      score_sig = if (all(is.na(r[2, ]))) 0 else max(r[2, ], na.rm = TRUE),
      score_cont = if (all(is.na(r[3, ]))) -1 else max(r[3, ], na.rm = TRUE),
      state_cor = sc)
  })
  cbind(pairs[, c("a", "b")], as.data.frame(do.call(rbind, res)))
}

## ---- metrics ---------------------------------------------------------------
auroc <- function(score, truth) {
  k <- !is.na(score); score <- score[k]; truth <- truth[k]   # failed pairs are NA -> dropped
  if (all(truth) || !any(truth)) return(NA_real_)
  r <- rank(score)
  (sum(r[truth]) - sum(truth) * (sum(truth) + 1) / 2) / (sum(truth) * sum(!truth))
}

## best TPR over thresholds whose FPR <= max_fpr
tpr_at_fpr <- function(score, truth, max_fpr = 0.05) {
  k <- !is.na(score); score <- score[k]; truth <- truth[k]
  if (all(truth) || !any(truth)) return(NA_real_)
  th <- sort(unique(score), decreasing = TRUE)
  ok <- vapply(th, function(t) mean(score[!truth] >= t) <= max_fpr, logical(1L))
  if (!any(ok)) return(0)
  max(vapply(th[ok], function(t) mean(score[truth] >= t), numeric(1L)))
}

metrics <- function(res, pairs, sname, arm) {
  tp <- pairs$truth
  m <- c(tpr = if (any(tp)) mean(res$called[tp] == 1, na.rm = TRUE) else NA,
         fpr = if (any(!tp)) mean(res$called[!tp] == 1, na.rm = TRUE) else NA,
         auroc_sig = auroc(res$score_sig, tp), auroc_cont = auroc(res$score_cont, tp),
         tpr05_sig = tpr_at_fpr(res$score_sig, tp), tpr05_cont = tpr_at_fpr(res$score_cont, tp),
         failed = mean(is.na(res$called)),
         n_false_shared = if (sname %in% all_private) sum(res$called, na.rm = TRUE) else NA,
         state_cor = if (any(tp)) mean(res$state_cor[tp], na.rm = TRUE) else NA)
  m[is.nan(m)] <- NA
  m
}

## ---- main ------------------------------------------------------------------
rows <- list()
add_rows <- function(sname, mf, arm, m) {
  m <- m[!is.na(m)]
  if (length(m) == 0L) return(invisible())
  rows[[length(rows) + 1L]] <<- data.frame(seed = seed, scenario = sname, missing_frac = mf,
                                           arm = arm, metric = names(m), value = unname(m))
}
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; list(v = v, s = proc.time()[["elapsed"]] - t0) }

for (sname in names(scen)) for (mf in c(0, 0.3)) {
  message(sprintf("[seed %d] %s missing=%.1f", seed, sname, mf))
  sim <- simulate_programs(scenario = scen[[sname]][1], N = N, G = G, missing_frac = mf,
                           seed = seed * 1000L + match(sname, names(scen)))
  nuis <- NULL
  if (scen[[sname]][2] == 1) {
    cf <- add_nuisance(sim$x, seed = seed * 1000L + 500L + match(sname, names(scen)))
    sim$x <- cf$x; nuis <- cf$nuis
  }
  ## the nuisance-fix arms run only on the confounded scenarios (S8, S9)
  tr <- list(z = lapply(sim$truth$activities, function(a) stats::setNames(a$z, paste0("obs", seq_along(a$z)))))
  pairs <- truth_pairs(sim$truth, names(sim$x$matrices))

  f1 <- timed(suppressWarnings(fit_celltype_programs(sim$x, loading_prior = "point_laplace",
                max_factors = 10L, var_type = 1L, backfit = TRUE, seed = seed, verbose = 0)))
  fit <- canonicalize_programs(f1$v)
  add_rows(sname, mf, "stage1", c(seconds = f1$s))

  arms <- list(
    PA_scores   = function() run_pa(fit, pairs, "scores", tr),
    PA_loadings = function() run_pa(fit, pairs, "loadings", tr),
    PMD_raw     = function() run_pmd(fit, pairs, tr))
  if (!is.null(nuis)) {
    fit_res <- residualize(fit, nuis)
    arms$PMD_resid <- function() run_pmd(fit_res, pairs, tr)
    arms$PA_scores_resid <- function() run_pa(fit_res, pairs, "scores", tr)
    est <- est_nuisance(fit)                     # no truth used
    nuis_est_cor <- abs(stats::cor(est, nuis[names(est)]))
    fit_est <- residualize(fit, est)
    arms$PMD_estnuis <- function() run_pmd(fit_est, pairs, tr)
    arms$PA_scores_estnuis <- function() run_pa(fit_est, pairs, "scores", tr)
  }
  for (arm in names(arms)) {
    r <- timed(tryCatch(arms[[arm]](), error = function(e) { message("  ", arm, " failed: ", conditionMessage(e)); NULL }))
    if (is.null(r$v)) { add_rows(sname, mf, arm, c(failed = 1)); next }
    ex <- if (grepl("estnuis", arm)) c(nuis_est_cor = nuis_est_cor) else NULL
    add_rows(sname, mf, arm, c(metrics(r$v, pairs, sname, arm), ex, seconds = r$s))
  }
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(do.call(rbind, rows), outfile, row.names = FALSE)  # incremental save
}
message("done: ", outfile)
