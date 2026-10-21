# Coordinated multicellular programs (CMPs)
#
# A CMP is a set of programs from >= 2 cell types whose donor-level scores
# co-vary, i.e. that track one shared latent donor activity. Gene overlap is
# irrelevant to the call: an IFN program in CD4 T cells and an inflammatory
# program in NK cells form a CMP if their donor scores co-vary, whatever their
# genes. Gene overlap, themes and ligand-receptor compatibility are separate
# annotations (annotate_cmp_lr(), annotate_cmp_themes(), cmp_gene_overlap()).
#
# Method, selected and ablated on cmp-bench-v1.1 and confirmed on held-out
# cmp-bench-v1.2 (benchmarks/cmp_detection/):
#   1. Rank-first scores: each program's donor scores are rank-transformed to
#      normal scores, then residualized on the blocking strata (e.g. site).
#      Point-Laplace priors produce single-donor spike factors that break
#      Pearson-based max-T calibration; ranks remove that failure mode.
#   2. Ambient filter: programs whose loadings track another cell type's
#      expression specificity (ambient RNA / doublet leakage) are dropped.
#   3. Edge statistic T = atanh(r) sqrt(n - 3) for every cross-cell-type
#      program pair over the donors the two cell types share.
#   4. Node statistic A_v = sum over partner cell types c of
#      max(0, max_{u in c} |T_vu| - tau)^2: a program's evidence is pooled
#      across partner cell types.
#   5. Null: donor rows permuted within strata, independently per cell type
#      (removes cross-cell-type coupling, keeps within-cell-type structure).
#      Program p-values are adjusted by the permutation maximum of A
#      (single-step max-T), so P(any CMP | no coordination) <= alpha.
#   6. Significant programs are linked by edges with |T| >= link_threshold
#      and grouped greedily, strongest edge first, one program per cell type.

.cmp_resolve_var <- function(fit, v, what, allow_na = FALSE) {
  if (is.null(v)) return(NULL)
  ids <- fit$observation_ids
  if (is.character(v) && length(v) == 1L && is.null(names(v))) {
    md <- fit$sample_metadata
    if (!v %in% names(md)) .stopf("%s column '%s' not found in sample_metadata.", what, v)
    out <- md[[v]]
    names(out) <- md$observation_id
  } else {
    if (is.null(names(v))) .stopf("%s must be a sample_metadata column name or a vector named by observation_id.", what)
    out <- v
  }
  miss <- setdiff(ids, names(out))
  if (length(miss)) .stopf("%s is missing for %d observation(s), e.g. '%s'.", what, length(miss), miss[1])
  out <- out[ids]
  if (!allow_na && anyNA(out)) .stopf("%s contains NA values.", what)
  out
}

# permutation / residualization blocks: strata, or strata x context in "within" mode
.cmp_blocks <- function(ids, strata, context, mode) {
  s <- if (is.null(strata)) rep("all", length(ids)) else as.character(strata[ids])
  if (mode == "within") {
    if (is.null(context)) .stopf("mode = 'within' requires `context`.")
    cx <- context[ids]
    lev <- if (is.numeric(cx) && length(unique(cx)) > 5L) {
      as.character(cut(cx, unique(stats::quantile(cx, 0:5 / 5)), include.lowest = TRUE))
    } else as.character(cx)
    s <- paste(s, lev, sep = ":")
  }
  stats::setNames(s, ids)
}

.cmp_rank_gauss <- function(S) {
  # NA scores stay NA: rank()'s default na.last=TRUE would assign them the last
  # (i.e. top) ranks and manufacture fake high-scoring donors (audit B-6).
  R <- apply(S, 2, function(v) {
    obs <- !is.na(v)
    r <- numeric(length(v))
    r[obs] <- stats::qnorm((rank(v[obs], ties.method = "average") - 0.5) / sum(obs))
    r[!obs] <- NA_real_
    r
  })
  if (is.null(dim(R))) R <- matrix(R, nrow(S))
  dimnames(R) <- dimnames(S)
  R
}

.cmp_prep <- function(scores, blocks, min_donors) {
  # Count OBSERVED donors, not rows: score matrices carry NA rows for donors
  # never observed in that cell type, so nrow() overstates the usable n.
  scores <- scores[vapply(scores, function(S) {
    !is.null(S) && ncol(S) > 0L && sum(rowSums(!is.na(S)) > 0L) >= min_donors
  }, TRUE)]
  lapply(scores, function(S) {
    R0 <- .cmp_rank_gauss(S)
    f <- factor(blocks[rownames(S)])
    X <- if (nlevels(f) > 1L) stats::model.matrix(~f) else matrix(1, nrow(S), 1L)
    # residualize each program on its observed donors; NA rows stay NA (audit B-6)
    R <- matrix(NA_real_, nrow(R0), ncol(R0), dimnames = dimnames(R0))
    for (j in seq_len(ncol(R0))) {
      obs <- !is.na(R0[, j])
      if (!any(obs)) next
      b <- qr.solve(X[obs, , drop = FALSE], R0[obs, j])
      R[obs, j] <- R0[obs, j] - X[obs, , drop = FALSE] %*% b
    }
    dimnames(R) <- dimnames(S)
    R
  })
}

.cmp_pairs <- function(Z, min_shared) {
  cts <- names(Z)
  out <- list()
  if (length(cts) < 2L) return(out)
  for (i in seq_len(length(cts) - 1L)) for (j in (i + 1L):length(cts)) {
    ids <- intersect(rownames(Z[[cts[i]]]), rownames(Z[[cts[j]]]))
    if (length(ids) >= min_shared) out[[length(out) + 1L]] <- list(a = cts[i], b = cts[j], ids = ids)
  }
  out
}

.cmp_edges <- function(Z, pairs) {
  do.call(rbind, lapply(pairs, function(p) {
    Za <- Z[[p$a]][p$ids, , drop = FALSE]
    Zb <- Z[[p$b]][p$ids, , drop = FALSE]
    # pairwise-complete donors per program pair; edge `n` reports the honest
    # pairwise-complete count rather than the nominal donor count (audit B-6)
    r <- suppressWarnings(stats::cor(Za, Zb, use = "pairwise.complete.obs"))
    n <- crossprod(!is.na(Za), !is.na(Zb))
    r[!is.finite(r)] <- 0
    t <- atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(pmax(n - 3, 0))
    t[n < 4] <- 0
    data.frame(ct_a = p$a, ct_b = p$b,
               prog_a = rownames(r)[row(r)], prog_b = colnames(r)[col(r)],
               stat = as.vector(abs(t)), sign = as.vector(sign(t)), n = as.vector(n),
               stringsAsFactors = FALSE)
  }))
}

.cmp_perm <- function(Z, blocks) {
  lapply(Z, function(S) {
    P <- S[.perm_within(nrow(S), blocks[rownames(S)]), , drop = FALSE]
    rownames(P) <- rownames(S)
    P
  })
}

.cmp_node_index <- function(edges) {
  ka <- paste(edges$ct_a, edges$prog_a, sep = "|")
  kb <- paste(edges$ct_b, edges$prog_b, sep = "|")
  grp <- factor(c(paste(ka, edges$ct_b, sep = "@"), paste(kb, edges$ct_a, sep = "@")))
  list(grp = grp, node = factor(sub("@[^@]*$", "", levels(grp))), ka = ka, kb = kb)
}

.cmp_node_stat <- function(stat, idx, tau) {
  m <- vapply(split(c(stat, stat), idx$grp), max, 1)
  v <- pmax(0, m - tau)^2
  out <- rowsum(v, idx$node, reorder = TRUE)[, 1]
  stats::setNames(as.numeric(out), names(out))
}

# greedy single linkage, strongest edge first, at most one program per cell type per group
.cmp_group <- function(edges) {
  edges <- edges[order(-edges$stat), , drop = FALSE]
  ka <- paste(edges$ct_a, edges$prog_a, sep = "|")
  kb <- paste(edges$ct_b, edges$prog_b, sep = "|")
  keys <- unique(c(ka, kb))
  cl <- stats::setNames(seq_along(keys), keys)
  for (i in seq_len(nrow(edges))) {
    ca <- cl[[ka[i]]]; cb <- cl[[kb[i]]]
    if (ca == cb) next
    ma <- names(cl)[cl == ca]; mb <- names(cl)[cl == cb]
    if (length(intersect(sub("\\|.*", "", ma), sub("\\|.*", "", mb)))) next
    cl[mb] <- ca
  }
  tab <- table(cl)
  keep <- as.integer(names(tab)[tab >= 2L])
  if (!length(keep)) return(NULL)
  cl <- cl[cl %in% keep]
  strength <- vapply(keep, function(k) max(edges$stat[ka %in% names(cl)[cl == k] & kb %in% names(cl)[cl == k]]), 1)
  rank_k <- keep[order(-strength)]
  ord <- order(match(cl, rank_k))
  cl <- cl[ord]
  data.frame(cmp = match(cl, rank_k), node = names(cl), stringsAsFactors = FALSE)
}

.cmp_empty <- function() {
  data.frame(cmp = integer(), cell_type = character(), program_id = character(), node_p = numeric(),
             sign = numeric(), ambient_flag = logical(), stringsAsFactors = FALSE)
}

#' Find coordinated multicellular programs
#'
#' Detects groups of programs from different cell types whose donor-level
#' scores co-vary (coordinated multicellular programs, CMPs), with
#' family-wise error control by within-strata donor permutation. Gene overlap
#' plays no role in the call; programs in a CMP may be thematically unrelated.
#'
#' @param fit A canonicalized `cell_program_fit` (see `canonicalize_programs()`).
#' @param strata Optional blocking variable (e.g. site or batch): a
#'   `sample_metadata` column name or a vector named by observation ID. Scores
#'   are residualized on it and donors are permuted only within it.
#' @param context Optional donor-level variable (e.g. disease status) used by
#'   `mode = "within"`; numeric variables with more than 5 distinct values are
#'   binned into quintiles.
#' @param mode `"marginal"` (default): any shared donor variation counts,
#'   including variation driven by `context`. `"within"`: only coordination that
#'   persists within levels of `context` (residualized on, and permuted within,
#'   strata x context).
#' @param ambient Treatment of ambient-like programs (see
#'   `ambient_program_scores()`): `"drop"` (default; excluded before testing),
#'   `"flag"` (tested, flagged in the output) or `"none"`.
#' @param ambient_threshold Ambient score above which a program is ambient-like
#'   (default 0.3).
#' @param profiles Optional per-cell-type mean log-expression profiles for the
#'   ambient score; defaults to column means of `fit$pseudobulk`.
#' @param alpha Family-wise error level (default 0.05).
#' @param n_perm Number of permutations (default 199).
#' @param tau Node-statistic offset; only edges with |T| > `tau` contribute
#'   (default 2).
#' @param link_threshold Minimum |T| for an edge to join two significant
#'   programs into one CMP (default 3).
#' @param min_shared Minimum donors shared by two cell types for their programs
#'   to be compared (default 20).
#' @param min_donors Minimum donors for a cell type to be used (default 20).
#' @param seed Random seed.
#'
#' @section Program cap:
#' CMP detection depends on each member program existing upstream. In the
#' cmp-bench-v1.2 held-out benchmark, fits with `max_factors = 30` raised
#' pooled membership F1 from 0.571 to 0.705 over the default `max_factors = 10`
#' (paired difference 0.127, 95% CI 0.113-0.142), with all false-call
#' guardrails <= 0.06. The `fit_celltype_programs()` default of 10 was tuned for
#' a different goal (donor-level disease signal in the concatenated
#' representation), so it is not changed; use `max_factors = 30` when the fit
#' is intended for `find_coordinated_programs()`. A one-time message is shown
#' when a fit with `max_factors < 30` has a cell type at its cap.
#'
#' @section Missing donor scores:
#' Programs with missing donor scores (NA) are handled by pairwise-complete
#' donors: rank transformation and strata residualization keep NA donors as NA,
#' cross-cell-type correlations use the donors observed in both programs, and
#' edge `n` reports the pairwise-complete donor count for that program pair. A
#' warning lists the affected programs. Missingness reduces power (fewer
#' effective donors); it does not manufacture high-scoring donors.
#'
#' @return An object of class `cmp_result`, a list with
#'   \describe{
#'     \item{cmps}{one row per (CMP, program): `cmp`, `cell_type`,
#'       `program_id`, `node_p` (FWER-adjusted), `sign` (orientation relative
#'       to the CMP's first program, propagated along linking edges as the
#'       product of edge signs along the path; CAND-05) and `ambient_flag`.}
#'     \item{edges}{edges among members of the same CMP (`stat` = |T|).}
#'     \item{nodes}{all tested programs with node statistic, adjusted p-value
#'       and ambient score.}
#'     \item{ambient}{ambient scores of all programs (or NULL).}
#'     \item{settings}{the arguments used.}
#'   }
#' @examples
#' # Toy canonicalized fit: program 1 of cell types A and B shares donor
#' # activity z (a planted coordinated program) and top genes G1-G15.
#' set.seed(1)
#' ids <- paste0("d", seq_len(40))
#' genes <- paste0("G", seq_len(100))
#' z <- rnorm(40)
#' mk <- function(ct) {
#'   S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
#'   S[, 1] <- z + rnorm(40, sd = 0.3)
#'   W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
#'               dimnames = list(genes, colnames(S)))
#'   W[1:15, 1] <- W[1:15, 1] + 1
#'   list(S = S, W = W)
#' }
#' ta <- mk("A"); tb <- mk("B")
#' fit <- structure(list(
#'   scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
#'   cell_types = c("A", "B"), observation_ids = ids,
#'   sample_metadata = data.frame(observation_id = ids,
#'                                group = rep(c("ctrl", "stim"), each = 20))),
#'   class = c("cell_program_fit", "cell_program_data", "list"))
#' res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
#' res$cmps
#' @export
find_coordinated_programs <- function(fit, strata = NULL, context = NULL,
                                      mode = c("marginal", "within"),
                                      ambient = c("drop", "flag", "none"),
                                      ambient_threshold = 0.3, profiles = NULL,
                                      alpha = 0.05, n_perm = 199L, tau = 2,
                                      link_threshold = 3, min_shared = 20L,
                                      min_donors = 20L, seed = 1L) {
  mode <- match.arg(mode)
  ambient <- match.arg(ambient)
  if (is.null(fit$scores) || is.null(fit$loadings)) {
    .stopf("`fit` has no scores/loadings; run canonicalize_programs() first.")
  }
  if ((n_perm + 1) * alpha < 1) .stopf("n_perm = %d is too small for alpha = %g.", n_perm, alpha)
  .cmp_cap_message(fit)
  ids <- fit$observation_ids
  st <- .cmp_resolve_var(fit, strata, "strata")
  cx <- .cmp_resolve_var(fit, context, "context")
  blocks <- .cmp_blocks(ids, st, cx, mode)
  settings <- list(strata = if (is.character(strata) && length(strata) == 1L) strata else if (is.null(strata)) NULL else "<vector>",
                   context = if (is.character(context) && length(context) == 1L) context else if (is.null(context)) NULL else "<vector>",
                   mode = mode, ambient = ambient, ambient_threshold = ambient_threshold, alpha = alpha,
                   n_perm = n_perm, tau = tau, link_threshold = link_threshold,
                   min_shared = min_shared, min_donors = min_donors, seed = seed)

  amb <- if (ambient != "none") ambient_program_scores(fit, profiles) else NULL
  Z <- .cmp_prep(fit$scores, blocks, min_donors)
  if (ambient == "drop" && !is.null(amb)) {
    bad <- amb[amb$ambient_score > ambient_threshold, , drop = FALSE]
    for (i in seq_len(nrow(bad))) {
      ct <- bad$cell_type[i]
      if (!is.null(Z[[ct]])) Z[[ct]] <- Z[[ct]][, colnames(Z[[ct]]) != bad$program_id[i], drop = FALSE]
    }
    Z <- Z[vapply(Z, ncol, 1L) > 0L]
  }
  naf <- do.call(rbind, lapply(names(Z), function(ct) {
    data.frame(cell_type = ct, program_id = colnames(Z[[ct]]), n_missing = colSums(is.na(Z[[ct]])),
               stringsAsFactors = FALSE)
  }))
  naf <- naf[naf$n_missing > 0L, , drop = FALSE]
  if (nrow(naf)) {
    warning(sprintf(paste0("find_coordinated_programs(): %d program(s) have missing donor scores (up to %d of %d donors, e.g. %s/%s); ",
                           "correlations use pairwise-complete donors and edge `n` reports the pairwise-complete count."),
                    nrow(naf), max(naf$n_missing), length(ids), naf$cell_type[1], naf$program_id[1]), call. = FALSE)
  }
  res <- structure(list(cmps = .cmp_empty(), edges = NULL, nodes = NULL, ambient = amb, settings = settings),
                   class = "cmp_result")
  pairs <- .cmp_pairs(Z, min_shared)
  if (!length(pairs)) return(res)

  withr::local_seed(seed)
  edges <- .cmp_edges(Z, pairs)
  idx <- .cmp_node_index(edges)
  A <- .cmp_node_stat(edges$stat, idx, tau)
  mx <- vapply(seq_len(n_perm), function(b) {
    max(.cmp_node_stat(.cmp_edges(.cmp_perm(Z, blocks), pairs)$stat, idx, tau))
  }, 1)
  pA <- vapply(A, function(a) (1 + sum(mx >= a)) / (n_perm + 1), 1)

  nodes <- data.frame(cell_type = sub("\\|.*", "", names(A)), program_id = sub("^[^|]*\\|", "", names(A)),
                      node_stat = unname(A), node_p = unname(pA), stringsAsFactors = FALSE)
  if (!is.null(amb)) {
    nodes$ambient_score <- amb$ambient_score[match(paste(nodes$cell_type, nodes$program_id),
                                                   paste(amb$cell_type, amb$program_id))]
  }
  res$nodes <- nodes
  sig <- names(pA)[pA < alpha]
  keep <- idx$ka %in% sig & idx$kb %in% sig & edges$stat >= link_threshold
  if (!any(keep)) return(res)
  g <- .cmp_group(edges[keep, , drop = FALSE])
  if (is.null(g)) return(res)

  cmps <- data.frame(cmp = g$cmp, cell_type = sub("\\|.*", "", g$node),
                     program_id = sub("^[^|]*\\|", "", g$node), node_p = unname(pA[g$node]),
                     stringsAsFactors = FALSE)
  key <- g$node
  cmps$sign <- NA_real_
  # CAND-05: orientation propagates along the LINKING edges (the kept edge set
  # that formed the CMP): a member's sign is the product of edge signs along
  # its path from the reference member. The previous lookup read the sign of
  # the DIRECT ref--member edge from the FULL edge table; for chain-joined
  # members that edge can be a non-linking noise edge, giving a random +/-1
  # orientation that propagated into cmp_scores/annotate_cmp_tf/cmp_triage.
  kn_a <- idx$ka[keep]; kn_b <- idx$kb[keep]; kn_s <- edges$sign[keep]
  for (k in unique(cmps$cmp)) {
    ii <- which(cmps$cmp == k)
    ref <- key[ii[1]]
    cmps$sign[ii[1]] <- 1
    oriented <- ref
    signs <- 1
    pending <- ii[-1]
    while (length(pending)) {
      hit <- vapply(pending, function(i) {
        m <- key[i]
        any((kn_a == m & kn_b %in% oriented) | (kn_b == m & kn_a %in% oriented))
      }, logical(1))
      if (!any(hit)) break  # defensive: .cmp_group components are connected
      i <- pending[which(hit)[1L]]
      m <- key[i]
      e <- which((kn_a == m & kn_b %in% oriented) |
                   (kn_b == m & kn_a %in% oriented))[1L]
      u <- if (kn_a[e] %in% oriented) kn_a[e] else kn_b[e]
      cmps$sign[i] <- signs[match(u, oriented)] * kn_s[e]
      oriented <- c(oriented, m)
      signs <- c(signs, cmps$sign[i])
      pending <- setdiff(pending, i)
    }
  }
  cmps$ambient_flag <- if (is.null(amb)) NA else
    amb$ambient_score[match(paste(cmps$cell_type, cmps$program_id), paste(amb$cell_type, amb$program_id))] > ambient_threshold
  ca <- cmps$cmp[match(idx$ka, key)]; cb <- cmps$cmp[match(idx$kb, key)]
  same <- !is.na(ca) & !is.na(cb) & ca == cb
  e_out <- edges[same, c("ct_a", "prog_a", "ct_b", "prog_b", "stat", "sign", "n"), drop = FALSE]
  e_out$cmp <- ca[same]
  rownames(cmps) <- NULL
  rownames(e_out) <- NULL
  res$cmps <- cmps
  res$edges <- e_out
  res
}

#' Print a CMP detection result
#'
#' @param x A `cmp_result` from [find_coordinated_programs()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.cmp_result <- function(x, ...) {
  ks <- unique(x$cmps$cmp)
  cat(sprintf("<cmp_result> %d coordinated multicellular program(s) at FWER %g (mode: %s, ambient: %s)\n",
              length(ks), x$settings$alpha, x$settings$mode, x$settings$ambient))
  for (k in ks) {
    m <- x$cmps[x$cmps$cmp == k, ]
    cat(sprintf("  CMP %d: %s\n", k, paste0(m$program_id, ifelse(!is.na(m$sign) & m$sign < 0, " (-)", ""),
                                            collapse = ", ")))
  }
  invisible(x)
}

#' Ambient-likeness of programs
#'
#' For each program, the largest absolute correlation between its gene loadings
#' and another cell type's expression specificity (`profile_other -
#' profile_own`, log scale). High values indicate a program that tracks
#' contamination by another cell type's transcripts (ambient RNA, doublets)
#' rather than biology intrinsic to the cell type. Genuine biology can also
#' score high (e.g. a monocyte-like program in DCs), so treat it as a flag.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param profiles Optional named list of per-cell-type mean log-expression
#'   vectors named by gene; defaults to column means of `fit$pseudobulk`.
#' @param min_genes Minimum shared genes to compare two cell types (default 50).
#' @return Data frame with `cell_type`, `program_id`, `ambient_score` and
#'   `source` (the cell type giving the maximum).
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
#' fit <- canonicalize_programs(fit)
#' ambient_program_scores(fit)
#' @keywords internal
ambient_program_scores <- function(fit, profiles = NULL, min_genes = 50L) {
  if (is.null(profiles)) {
    if (is.null(fit$pseudobulk)) .stopf("No `profiles` given and `fit$pseudobulk` is absent.")
    profiles <- lapply(fit$pseudobulk, function(Y) colMeans(as.matrix(Y), na.rm = TRUE))
  }
  out <- lapply(fit$cell_types, function(ct) {
    W <- fit$loadings[[ct]]
    if (is.null(W) || !ncol(W) || is.null(profiles[[ct]])) return(NULL)
    best <- rep(0, ncol(W)); src <- rep(NA_character_, ncol(W))
    for (s in setdiff(names(profiles), ct)) {
      g <- intersect(rownames(W), intersect(names(profiles[[s]]), names(profiles[[ct]])))
      if (length(g) < min_genes) next
      spec <- profiles[[s]][g] - profiles[[ct]][g]
      r <- abs(suppressWarnings(stats::cor(W[g, , drop = FALSE], spec)))[, 1]
      r[!is.finite(r)] <- 0
      up <- r > best
      best[up] <- r[up]; src[up] <- s
    }
    data.frame(cell_type = ct, program_id = colnames(W), ambient_score = best, source = src,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(cell_type = character(), program_id = character(), ambient_score = numeric(),
                      source = character(), stringsAsFactors = FALSE)
  }
  out
}

#' Consensus donor scores of coordinated programs
#'
#' Mean of the member programs' standardized scores, each oriented by its sign
#' within the CMP (available-case over cell types for each donor).
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` data frame.
#' @return Matrix, donors x CMPs.
#' @examples
#' # Toy canonicalized fit: program 1 of cell types A and B shares donor
#' # activity z (a planted coordinated program) and top genes G1-G15.
#' set.seed(1)
#' ids <- paste0("d", seq_len(40))
#' genes <- paste0("G", seq_len(100))
#' z <- rnorm(40)
#' mk <- function(ct) {
#'   S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
#'   S[, 1] <- z + rnorm(40, sd = 0.3)
#'   W <- matrix(rnorm(length(genes) * 2, sd = 0.1), length(genes), 2,
#'               dimnames = list(genes, colnames(S)))
#'   W[1:15, 1] <- W[1:15, 1] + 1
#'   list(S = S, W = W)
#' }
#' ta <- mk("A"); tb <- mk("B")
#' fit <- structure(list(
#'   scores = list(A = ta$S, B = tb$S), loadings = list(A = ta$W, B = tb$W),
#'   cell_types = c("A", "B"), observation_ids = ids,
#'   sample_metadata = data.frame(observation_id = ids,
#'                                group = rep(c("ctrl", "stim"), each = 20))),
#'   class = c("cell_program_fit", "cell_program_data", "list"))
#' res <- find_coordinated_programs(fit, n_perm = 99, ambient = "none")
#' S <- cmp_scores(fit, res)
#' head(S)
#' @keywords internal
cmp_scores <- function(fit, cmp) {
  tab <- if (inherits(cmp, "cmp_result")) cmp$cmps else cmp
  ks <- unique(tab$cmp)
  ids <- fit$observation_ids
  M <- matrix(NA_real_, length(ids), length(ks), dimnames = list(ids, paste0("CMP", ks)))
  for (j in seq_along(ks)) {
    m <- tab[tab$cmp == ks[j], , drop = FALSE]
    V <- matrix(NA_real_, length(ids), nrow(m), dimnames = list(ids, NULL))
    for (i in seq_len(nrow(m))) {
      v <- fit$scores[[m$cell_type[i]]][, m$program_id[i]]
      sg <- if (!is.null(m$sign) && is.finite(m$sign[i])) m$sign[i] else 1
      V[names(v), i] <- sg * as.numeric(scale(v))
    }
    M[, j] <- rowMeans(V, na.rm = TRUE)
    M[rowSums(!is.na(V)) == 0, j] <- NA
  }
  M
}

## cmp-bench-v1.2: recall of CMP member programs saturates at max_factors = 30; the fitting default (10)
## was tuned for disease-signal recovery. Message once per session, only if some cell type hit the cap.
.cmp_cap_message <- function(fit) {
  k <- fit$preprocessing$max_factors
  if (!is.numeric(k) || length(k) != 1L || k >= 30 || isTRUE(getOption("cellprograms.cmp_cap_message_shown"))) return(invisible(FALSE))
  at_cap <- names(Filter(function(S) !is.null(S) && ncol(S) >= k, fit$scores))
  if (!length(at_cap)) return(invisible(FALSE))
  options(cellprograms.cmp_cap_message_shown = TRUE)
  message(sprintf(paste0("find_coordinated_programs(): %d cell type(s) reached max_factors = %g (e.g. %s). ",
                         "For CMP detection, refitting with max_factors = 30 recovered more member programs in the ",
                         "cmp-bench-v1.2 benchmark. Shown once per session."), length(at_cap), k, at_cap[1]))
  invisible(TRUE)
}
