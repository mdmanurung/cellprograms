## score.R — scorer for the coordination benchmark.
##
## Contract: every method returns clusters = data.frame(cluster, cell_type,
## program_id, ...) and optionally activity (donors x clusters). A
## "coordination call" is a cluster spanning >= 2 distinct cell types.
##
## Matching: one-to-one Hungarian assignment (clue::solve_LSAP) maximizing
## total membership-set F1 between calls and planted recurring programs.
## An assignment is ACCEPTED iff recall >= 0.5 AND precision >= 0.5 (the
## prereg "found" rule); rejected/unassigned calls are spurious. Spurious
## calls with donor activity correlated with the nuisance (|r| >= 0.5) are
## decoy calls. In priv_corr mode a spurious call containing both members of
## a correlated-private pair is a trap failure.

.f1 <- function(a, b) {
  inter <- length(intersect(a, b))
  if (inter == 0L) return(0)
  2 * inter / (length(a) + length(b))
}

score_method <- function(res, truth, decoy_cor = 0.5) {
  cl <- res$clusters
  progs <- truth$programs
  call_sets <- if (nrow(cl)) split(cl$cell_type, cl$cluster) else list()
  call_sets <- lapply(call_sets, unique)
  call_sets <- call_sets[lengths(call_sets) >= 2L]     # coordination calls only
  n_calls <- length(call_sets)

  prog_level <- data.frame()
  accepted_call <- rep(FALSE, n_calls)
  if (length(progs)) {
    P <- length(progs)
    F1 <- matrix(0, n_calls, P)
    for (i in seq_len(n_calls)) for (j in seq_len(P)) {
      F1[i, j] <- .f1(call_sets[[i]], progs[[j]]$members)
    }
    ## pad to square for rectangular cases; dummies have F1 = 0
    n <- max(n_calls, P, 1L)
    F1p <- matrix(0, n, n)
    if (n_calls > 0L) F1p[seq_len(n_calls), seq_len(P)] <- F1
    asg <- if (n_calls > 0L) as.integer(clue::solve_LSAP(F1p, maximum = TRUE)) else rep(seq_len(n), 0)
    rows <- list()
    for (j in seq_len(P)) {
      i <- if (length(asg)) which(asg == j) else integer(0)
      i <- i[i <= n_calls]
      rec <- prec <- 0; be <- NA_integer_; found <- FALSE; cid <- NA_integer_
      if (length(i) == 1L) {
        cs <- call_sets[[i]]; tm <- progs[[j]]$members
        inter <- length(intersect(cs, tm))
        rec <- inter / length(tm); prec <- inter / length(cs)
        if (rec >= 0.5 && prec >= 0.5) {
          found <- TRUE; accepted_call[i] <- TRUE
          be <- length(cs) - length(tm); cid <- as.integer(names(call_sets)[i])
        }
      }
      rows[[j]] <- data.frame(program_id = progs[[j]]$id, n_true = length(progs[[j]]$members),
                              matched_cluster = cid, recall = rec, precision = prec,
                              breadth_error = be, found = found, stringsAsFactors = FALSE)
    }
    prog_level <- do.call(rbind, rows)
  }

  ## spurious / decoy / trap among calls
  spurious <- if (n_calls) which(!accepted_call) else integer(0)
  n_decoy <- 0L; n_trap <- 0L
  if (length(spurious)) {
    bz <- truth$nuisance
    for (i in spurious) {
      cid <- names(call_sets)[i]
      if (!is.null(res$activity) && paste0("cluster_", cid) %in% colnames(res$activity)) {
        a <- res$activity[, paste0("cluster_", cid)]
        d <- intersect(names(a)[is.finite(a)], names(bz))
        if (length(d) >= 20L && stats::sd(a[d]) > 0) {
          if (abs(stats::cor(a[d], bz[d])) >= decoy_cor) n_decoy <- n_decoy + 1L
        }
      }
      for (dp in truth$decoy_pairs) {
        if (all(dp$members %in% call_sets[[i]])) { n_trap <- n_trap + 1L; break }
      }
    }
  }
  dataset_level <- data.frame(
    method = res$method, n_calls = n_calls, n_spurious = length(spurious),
    n_decoy = n_decoy, n_trap = n_trap,
    spurious_per_call = if (n_calls) length(spurious) / n_calls else 0,
    stringsAsFactors = FALSE)
  list(program_level = prog_level, dataset_level = dataset_level)
}

## Self-check: hand-built cases (run by test_score.R)
score_selfcheck <- function() {
  tr <- list(programs = list(
    list(id = "R50", members = c("ct01", "ct03", "ct05", "ct08", "ct10")),
    list(id = "R20", members = c("ct02", "ct09"))),
    decoy_pairs = list(), nuisance = setNames(rnorm(100), sprintf("d%03d", 1:100)))
  ## perfect R50, partial R20 (1/2 members -> recall 0.5 prec 1 -> found), one spurious
  cl <- data.frame(cluster = c(1, 1, 1, 1, 1, 2, 2, 3, 3),
                   cell_type = c("ct01", "ct03", "ct05", "ct08", "ct10",
                                 "ct02", "ct09", "ct04", "ct06"),
                   program_id = NA, stat = 1, fwer_p = NA)
  s <- score_method(list(method = "chk", clusters = cl, activity = NULL), tr)
  stopifnot(s$program_level$found[1], s$program_level$found[2],
            s$dataset_level$n_calls == 3L, s$dataset_level$n_spurious == 1L)
  ## R20 with only 1 of 2 members in a 3-ct cluster: recall 0.5, prec 1/3 -> not found, spurious
  cl2 <- data.frame(cluster = c(1, 1, 1), cell_type = c("ct02", "ct04", "ct06"),
                    program_id = NA, stat = 1, fwer_p = NA)
  s2 <- score_method(list(method = "chk", clusters = cl2, activity = NULL), tr)
  stopifnot(!s2$program_level$found[s2$program_level$program_id == "R20"],
            s2$dataset_level$n_spurious == 1L)
  cat("score self-checks passed\n")
}
