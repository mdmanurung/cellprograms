# program-sharing.R — post-hoc shared-vs-specific program assessment (mashr)
#
# Doctrine: programs are fitted independently per cell type; cross-cell-type
# sharing of GENE LOADINGS is assessed post-hoc, never assumed. Shared loadings
# are a different axis from coordinated DONOR ACTIVITY: two cell types can use
# the same biological program (shared gene loadings) with uncorrelated donor
# scores, and vice versa. Donor-score coordination (CMP detection) is a
# separate workstream — see find_coordinated_programs(); program_sharing()
# says nothing about donor-score coordination.
#
# mashr mapping (Urbut et al., Nat Genet 2019): units = genes on the shared
# panel; conditions = cell types; effects = canonicalized program loadings;
# SEs = loading posterior SDs (flashier F_psd, scale-canonicalized like the
# loadings; sign flips leave SDs unchanged). mashr is a Suggests dependency.

.has_mashr <- function() requireNamespace("mashr", quietly = TRUE)

# Validate a canonicalized fit; with psd = TRUE also require loading posterior
# SDs (old fits predate F_psd storage -> graceful refit error).
.sharing_check_fit <- function(fit, psd = TRUE) {
  if (!inherits(fit, "cell_program_fit")) .stopf("`fit` must be a cell_program_fit.")
  if (is.null(fit$loadings) || is.null(fit$canonicalization)) {
    .stopf("`fit` has no canonicalized loadings; run canonicalize_programs() first.")
  }
  if (psd) {
    for (ct in names(fit$loadings)) {
      if (ncol(fit$loadings[[ct]]) == 0L) next
      f <- fit$fits[[ct]]$flash
      if (is.null(f$F_psd) || nrow(f$F_psd) != nrow(fit$loadings[[ct]])) {
        .stopf(paste0("Cell type '%s': loading posterior SDs (F_psd) are not stored in this fit. ",
                      "Refit with fit_celltype_programs() from the current cellprograms version ",
                      "(posterior-SD storage was added together with program_sharing)."), ct)
      }
    }
  }
  invisible(fit)
}

#' Numerically degenerate programs, two kinds:
#' 1. flashier can retain factors whose loadings underflow to numerical noise
#'    (median |loading| many orders of magnitude below the cell type's strongest
#'    program; e.g. sparse factors whose support is a few genes and whose
#'    remaining entries sit at ~1e-16). Their posterior SDs underflow as well,
#'    producing garbage z-scores, and correlation-based matching is
#'    scale-invariant, so they spuriously match every strong program.
#' 2. Fixed, never-updated factors: real-scale loadings whose posterior SD
#'    column collapsed to machine precision (observed for the smallest cell
#'    type in both real cohorts). Point estimates without usable uncertainty:
#'    mashr cannot assess them.
#' A program is degenerate when median |loading| < `rel_tol` x the cell type's
#' largest, or when its median posterior SD (scaled to the loadings) is <=
#' `sd_floor`.
#' @keywords internal
.degenerate_programs <- function(fit, rel_tol = 1e-6, sd_floor = 1e-12) {
  out <- list()
  for (ct in names(fit$loadings)) {
    W <- fit$loadings[[ct]]
    if (!ncol(W)) next
    sc <- apply(abs(W), 2L, stats::median)
    bad <- sc < rel_tol * max(sc)
    if (!is.null(fit$fits[[ct]]$flash$F_psd)) {
      sdm <- apply(fit$fits[[ct]]$flash$F_psd, 2L, stats::median) *
        fit$canonicalization[[ct]]$scale_multiplier
      bad <- bad | (sdm <= sd_floor)
    }
    if (any(bad)) out[[ct]] <- colnames(W)[bad]
  }
  out
}


#' Candidate cross-cell-type program matches from loading correlation
#'
#' For every pair of cell types, correlates (absolute Pearson) each program's
#' canonicalized loading vector over the shared gene universe (intersection of
#' the two feature sets) and calibrates the statistic against a permutation
#' null: gene labels are permuted within each cell type's loading matrix
#' (`n_perm` permutations; per-pair p = (1 + exceedances) / (n_perm + 1)).
#' `q` is BH across all tested program pairs.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param n_perm Permutations for the null (default 199; >= 1 required).
#' @param seed Master seed; each cell-type pair gets a deterministic sub-seed.
#' @return A data.frame with columns ct_a, prog_a, ct_b, prog_b, cor, p, q,
#'   sorted by q then -|cor|. Rows with q < alpha are candidate matches.
.match_programs_loadings <- function(fit, n_perm = 199L, seed = 1L) {
  .sharing_check_fit(fit, psd = FALSE)
  if (n_perm < 1L) .stopf("`n_perm` must be >= 1.")
  cts <- names(fit$loadings)
  cts <- cts[vapply(fit$loadings, function(W) ncol(W) > 0L, logical(1L))]
  if (length(cts) < 2L) .stopf("Need >= 2 cell types that retained a program; found %d.", length(cts))
  degen <- .degenerate_programs(fit)
  if (length(degen)) {
    .warnf("Excluding %d numerically degenerate program(s) (%s) from matching.",
           sum(lengths(degen)), paste(unlist(degen), collapse = ", "))
  }
  rows <- list()
  for (pr in utils::combn(cts, 2L, simplify = FALSE)) {
    Wa <- fit$loadings[[pr[1L]]]; Wb <- fit$loadings[[pr[2L]]]
    if (pr[1L] %in% names(degen)) Wa <- Wa[, setdiff(colnames(Wa), degen[[pr[1L]]]), drop = FALSE]
    if (pr[2L] %in% names(degen)) Wb <- Wb[, setdiff(colnames(Wb), degen[[pr[2L]]]), drop = FALSE]
    if (!ncol(Wa) || !ncol(Wb)) next
    genes <- intersect(rownames(Wa), rownames(Wb))
    if (length(genes) < 10L) {
      .warnf("Gene intersection for %s vs %s is %d (< 10); pair skipped.", pr[1L], pr[2L], length(genes))
      next
    }
    Wa <- Wa[genes, , drop = FALSE]; Wb <- Wb[genes, , drop = FALSE]
    obs <- stats::cor(Wa, Wb)
    exc <- matrix(0, nrow(obs), ncol(obs))
    set.seed(.pair_seed(seed, pr))
    G <- length(genes)
    for (p in seq_len(n_perm)) {
      null <- stats::cor(Wa[sample.int(G), , drop = FALSE],
                         Wb[sample.int(G), , drop = FALSE])
      exc <- exc + (abs(null) >= abs(obs))
    }
    pv <- (1 + exc) / (n_perm + 1)
    for (i in seq_len(ncol(Wa))) {
      for (j in seq_len(ncol(Wb))) {
        rows[[length(rows) + 1L]] <- data.frame(
          ct_a = pr[1L], prog_a = colnames(Wa)[i],
          ct_b = pr[2L], prog_b = colnames(Wb)[j],
          cor = obs[i, j], p = pv[i, j], stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) .stopf("No cell-type pair had a usable gene intersection (>= 10 genes).")
  out <- do.call(rbind, rows)
  out$q <- stats::p.adjust(out$p, "BH")
  out <- out[order(out$q, -abs(out$cor)), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_perm") <- n_perm
  out
}

# Greedy agglomeration of significant pairwise matches into disjoint matched
# groups with at most one program per cell type (edges in q / -|cor| order;
# an edge that would put two programs of the same cell type in one group is
# skipped). `matches` needs ct_a/prog_a/ct_b/prog_b; q optional (absent = all
# rows are edges).
.sharing_groups <- function(matches, alpha = 0.05) {
  if (!all(c("ct_a", "prog_a", "ct_b", "prog_b") %in% names(matches))) {
    .stopf("`matches` must have columns ct_a, prog_a, ct_b, prog_b.")
  }
  edges <- if ("q" %in% names(matches)) {
    matches[!is.na(matches$q) & matches$q < alpha, , drop = FALSE]
  } else {
    matches
  }
  if (!nrow(edges)) return(list())
  edges <- if ("q" %in% names(edges)) {
    edges[order(edges$q, -abs(edges$cor)), , drop = FALSE]
  } else {
    edges[order(-abs(edges$cor)), , drop = FALSE]
  }
  node_ct <- character(0); lab <- integer(0); members <- list()
  for (e in seq_len(nrow(edges))) {
    u <- paste0(edges$ct_a[e], "\r", edges$prog_a[e])
    v <- paste0(edges$ct_b[e], "\r", edges$prog_b[e])
    node_ct[u] <- edges$ct_a[e]; node_ct[v] <- edges$ct_b[e]
    lu <- if (u %in% names(lab)) lab[[u]] else NA_integer_
    lv <- if (v %in% names(lab)) lab[[v]] else NA_integer_
    if (is.na(lu) && is.na(lv)) {
      members[[length(members) + 1L]] <- c(u, v)
      lab[u] <- lab[v] <- length(members)
    } else if (!is.na(lu) && is.na(lv)) {
      if (edges$ct_b[e] %in% node_ct[members[[lu]]]) next
      members[[lu]] <- c(members[[lu]], v); lab[v] <- lu
    } else if (is.na(lu) && !is.na(lv)) {
      if (edges$ct_a[e] %in% node_ct[members[[lv]]]) next
      members[[lv]] <- c(members[[lv]], u); lab[u] <- lv
    } else if (lu != lv) {
      if (any(node_ct[members[[lu]]] %in% node_ct[members[[lv]]])) next
      members[[lu]] <- c(members[[lu]], members[[lv]])
      lab[members[[lv]]] <- lu
      members[[lv]] <- character(0)          # keep list indices stable
    }
  }
  groups <- lapply(members, function(nds) {
    nds <- unique(nds)
    if (length(nds) < 2L) return(NULL)
    data.frame(cell_type = unname(node_ct[nds]),
               program_id = sub("^[^\r]+\r", "", nds),
               stringsAsFactors = FALSE)
  })
  groups <- Filter(Negate(is.null), groups)
  names(groups) <- paste0("group", seq_along(groups))
  groups
}

# Coerce a matched-set specification to a data.frame(cell_type, program_id).
.as_group <- function(matches) {
  if (is.data.frame(matches) && all(c("cell_type", "program_id") %in% names(matches))) {
    g <- matches[, c("cell_type", "program_id")]
  } else if (is.data.frame(matches) &&
             all(c("ct_a", "prog_a", "ct_b", "prog_b") %in% names(matches)) && nrow(matches) == 1L) {
    g <- data.frame(cell_type = c(matches$ct_a, matches$ct_b),
                    program_id = c(matches$prog_a, matches$prog_b),
                    stringsAsFactors = FALSE)
  } else {
    .stopf("A matched set must be a data.frame with columns cell_type/program_id (one row per cell type), or a single ct_a/prog_a/ct_b/prog_b row.")
  }
  if (nrow(g) < 2L || anyDuplicated(g$cell_type)) {
    .stopf("A matched set needs >= 2 rows with distinct cell types.")
  }
  rownames(g) <- NULL
  g
}

# D-2a/D-2c: validate that every referenced (cell_type, program_id) exists in
# the fit BEFORE the mashr layer indexes into it — unknown references used to
# crash with opaque subscript errors deep in .sharing_matrices.
.validate_groups <- function(fit, groups) {
  tab <- .program_table(fit)
  have <- paste(tab$cell_type, tab$program_id, sep = "\r")
  for (gi in seq_along(groups)) {
    g <- groups[[gi]]
    want <- paste(g$cell_type, g$program_id, sep = "\r")
    bad <- !want %in% have
    if (any(bad)) {
      .stopf("matches[[%d]] references unknown program(s): %s.",
             gi, paste(sprintf("%s/%s", g$cell_type[bad], g$program_id[bad]), collapse = ", "))
    }
  }
  invisible(TRUE)
}

# G x C sign-aligned loading (B) and SE (S) matrices for one matched set over
# the shared gene panel. Loadings are already canonicalized within each cell
# type; members matched on |cor| are additionally sign-aligned to the first
# (reference) member. Posterior SDs get the canonicalization scale multiplier
# (SDs scale like loadings; sign flips leave them unchanged).
.sharing_matrices <- function(fit, group, min_genes = 50L) {
  .sharing_check_fit(fit, psd = TRUE)
  degen <- .degenerate_programs(fit)
  drop_j <- which(vapply(seq_len(nrow(group)), function(j) {
    ids <- degen[[group$cell_type[j]]]
    !is.null(ids) && group$program_id[j] %in% ids
  }, logical(1L)))
  if (length(drop_j)) {
    .warnf("Dropping %d matched condition(s) with numerically degenerate programs (%s).",
           length(drop_j), paste(group$program_id[drop_j], collapse = ", "))
    group <- group[-drop_j, , drop = FALSE]
    if (nrow(group) < 2L) {
      .stopf("Matched group has fewer than 2 non-degenerate programs after dropping degenerate members.")
    }
  }
  genes <- Reduce(intersect, lapply(group$cell_type, function(ct) rownames(fit$loadings[[ct]])))
  if (length(genes) < min_genes) {
    .stopf("Shared gene panel for {%s} is %d genes (< %d).",
           paste(group$cell_type, collapse = ", "), length(genes), min_genes)
  }
  C <- nrow(group)
  B <- S <- matrix(NA_real_, length(genes), C, dimnames = list(genes, group$cell_type))
  for (j in seq_len(C)) {
    ct <- group$cell_type[j]; prog <- group$program_id[j]
    W <- fit$loadings[[ct]]
    k <- match(prog, colnames(W))
    if (is.na(k)) .stopf("Unknown program '%s' in cell type '%s'.", prog, ct)
    s_k <- fit$canonicalization[[ct]]$scale_multiplier[k]
    B[, j] <- W[genes, k]
    S[, j] <- fit$fits[[ct]]$flash$F_psd[match(genes, rownames(W)), k] * s_k
  }
  # Degenerate posteriors from the sparse point-laplace prior: genes unused by a member's
  # program get loadings at the numerical noise floor (1e-16 and below, with an exponential
  # tail) AND underflowed SDs, which mashr's zero check rejects. Per mashr's own guidance,
  # set a small positive SE there: z ~ 0 -> lfsr = 1 -> the gene is (correctly) called
  # inactive in that condition and the panel is preserved. The floor is relative to the
  # condition's own loading scale (the SD scale itself can be collapsed by underflow).
  # Entries with a zero SD but a MATERIALLY nonzero loading (relative to that same scale;
  # an inconsistent posterior) or non-finite values are genuinely broken: drop.
  tiny <- S <= 1e-12
  b_scale <- pmax(apply(abs(B), 2L, stats::median), 1e-6)
  broken <- rowSums(!is.finite(B) | !is.finite(S) |
                      (tiny & abs(B) > 1e-6 * matrix(b_scale, nrow(B), ncol(B), byrow = TRUE))) > 0L
  if (any(broken)) {
    .warnf("Dropping %d/%d shared-panel genes with non-finite or inconsistent loadings/SDs.",
           sum(broken), length(genes))
    B <- B[!broken, , drop = FALSE]; S <- S[!broken, , drop = FALSE]
    genes <- genes[!broken]; tiny <- tiny[!broken, , drop = FALSE]
  }
  if (any(tiny)) {
    s_floor <- b_scale / 100
    S[tiny] <- matrix(rep(s_floor, each = nrow(S)), nrow(S), ncol(S))[tiny]
    .warnf(paste0("Set floor SEs on %d zero-SD loading entries (point-laplace noise-floor ",
                  "posteriors): mashr calls them inactive in those conditions."),
           sum(tiny))
  }
  if (nrow(B) < min_genes) {
    .stopf("Shared gene panel for {%s} has %d usable genes (< %d) after dropping degenerate entries.",
           paste(group$cell_type, collapse = ", "), nrow(B), min_genes)
  }
  for (j in seq_len(C)[-1L]) {
    r <- stats::cor(B[, 1L], B[, j])
    if (is.finite(r) && r < 0) B[, j] <- -B[, j]
  }
  list(B = B, S = S, genes = genes, cell_types = group$cell_type,
       programs = group$program_id)
}

#' CorShrink-style null-correlation estimate for a matched program set
#'
#' Estimates the correlation of loading estimation noise across cell types
#' (shared donors, ambient leakage) from the z-scores (loading / posterior SD)
#' of null genes — genes with |z| < `z_thresh` in EVERY condition of the
#' matched set — and shrinks the correlation matrix toward the identity with a
#' Ledoit-Wolf-style analytic intensity: lambda = sum of plug-in sampling
#' variances (1 - r^2)^2 / (n - 1) over the sum of squared off-diagonal
#' correlations, clipped to [0, 1]. With fewer than `min_null` null genes the
#' identity is returned (conservative; warns).
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param matches One matched set: a data.frame with columns `cell_type` and
#'   `program_id` (one row per cell type), or a single-row pairwise match with
#'   columns ct_a, prog_a, ct_b, prog_b.
#' @param z_thresh |z| below this in every condition marks a null gene.
#' @param min_null Minimum number of null genes for a usable estimate.
#' @return The shrunk C x C correlation matrix (dimnames = cell types), with
#'   attributes `n_null` and `lambda`.
#' @keywords internal
.estimate_V <- function(fit, matches, z_thresh = 1, min_null = 20L) {
  group <- .as_group(matches)
  mats <- .sharing_matrices(fit, group, min_genes = 10L)
  Z <- mats$B / mats$S
  C <- ncol(Z)
  null_g <- rowSums(abs(Z) >= z_thresh) == 0L
  lambda <- NA_real_
  if (sum(null_g) < min_null) {
    .warnf("Only %d null genes (|z| < %g in all %d conditions); using V = identity.",
           sum(null_g), z_thresh, C)
    V <- diag(C)
  } else {
    R <- stats::cor(Z[null_g, , drop = FALSE])
    if (anyNA(R)) {
      .warnf("Degenerate null z-scores; using V = identity.")
      R <- diag(C)
    }
    n <- sum(null_g)
    off <- R[upper.tri(R)]
    v <- (1 - off^2)^2 / (n - 1)
    lambda <- if (sum(off^2) == 0) 1 else min(1, sum(v) / sum(off^2))
    V <- (1 - lambda) * R + lambda * diag(C)
  }
  dimnames(V) <- list(group$cell_type, group$cell_type)
  attr(V, "n_null") <- sum(null_g)
  attr(V, "lambda") <- lambda
  V
}

#' Post-hoc shared-vs-specific assessment of fitted programs (mashr)
#'
#' Assesses which fitted programs are the SAME biological program across cell
#' types (shared gene loadings) versus cell-type-private, using mashr-style
#' adaptive shrinkage: for each matched program group (one program per cell
#' type) the G x C matrix of canonicalized loadings and their posterior SDs go
#' to [mashr::mash()] with the canonical plus data-driven (PCA/ED) covariance
#' set and the CorShrink-style null-correlation matrix [`.estimate_V`]. The
#' estimated mixture weights over covariance patterns ARE the empirical sharing
#' patterns (e.g. `equal_effects` = shared across all conditions,
#' `singletons_*` = one cell type only).
#'
#' This concerns shared GENE LOADINGS only. Coordinated DONOR ACTIVITY (CMP
#' detection, [find_coordinated_programs()]) is a separate axis and a separate
#' workstream; a shared program need not have correlated donor scores.
#'
#' @param fit A canonicalized `cell_program_fit` with stored loading posterior
#'   SDs (fits from older versions must be refit).
#' @param matches `NULL` (default): run [`.match_programs_loadings`] and group
#'   matches significant at FDR `alpha`. Alternatively a pairwise match table
#'   (columns ct_a/prog_a/ct_b/prog_b, optionally q), one matched set, or a
#'   list of matched sets (data.frames with columns cell_type, program_id).
#' @param alpha FDR level for candidate matches and the strong-gene subset for
#'   data-driven covariances (default 0.05).
#' @param n_top Number of top-|loading| genes per group for the summary
#'   fractions (default 50). Choose relative to the smallest shared support
#'   you want to detect: a group sharing S genes can reach
#'   `frac_top_shared` at most S/`n_top`, so with the 0.5 "shared" threshold
#'   any support smaller than `n_top`/2 (25 at the default) classifies at
#'   best "partial" — never "shared" — no matter how clean the signal
#'   (confirm-scale evidence: planted 10-gene sex/genotype-shared programs
#'   were never promoted in any arm; CP-4 Amendment 1, limitation register
#'   L3). For narrow shared modules (VDJ/sex/genotype scale, ~10-30 genes),
#'   lower `n_top` accordingly, and report `frac_top_shared`,
#'   `frac_top_private` and `p_all_shared` directly rather than the ternary
#'   label alone.
#' @param lfsr_thresh lfsr threshold calling a gene active in a condition
#'   (default 0.05).
#' @param n_perm,seed Permutation calibration for the matching step. CAND-07
#'   caveat: the gene-permutation null treats genes as exchangeable;
#'   block-correlated (co-regulated) genes violate that assumption, so
#'   matching q-values are anticonservative when a large correlated gene block
#'   drives both programs (audit D-1: 67.0% false-match rate under block
#'   equicorrelation 0.55 vs 2.0% iid). A block-aware null is pending a
#'   registered calibration benchmark; interpret marginal matches accordingly.
#' @param min_genes Minimum shared-panel size per group (default 50).
#' @param min_null Minimum null genes for V estimation (else identity).
#' @param z_thresh |z| below this in every condition marks a null gene for V
#'   estimation (default 1). With many conditions the intersection at 1 can be
#'   empty; raise to ~1.5 for large matched groups.
#' @param data_driven Add PCA/ED data-driven covariances to the canonical set
#'   (default TRUE; skipped with a warning when too few strong genes).
#' @param npc Number of PCA covariance components (capped at #conditions).
#' @param verbose Let mashr print progress (default FALSE).
#' @param ... Passed to [mashr::mash()].
#' @return An object of class `program_sharing`:
#'   * `matches`: the pairwise candidate-match table (NULL for explicit groups);
#'   * `groups`: list of matched-set data.frames (cell_type, program_id);
#'   * `results`: per-group list with `B_hat`/`S_hat` (mash inputs), `V`,
#'     `posterior_mean` and `lfsr` (G x C), `mixture_weights` (mashr's pi over
#'     the scaled covariance grid), `pattern_weights` (pi aggregated over base
#'     covariance patterns), `p_all_shared` (mean posterior weight of the
#'     equal_effects pattern over the top genes), `frac_top_shared` /
#'     `frac_top_private` (fraction of the top `n_top` genes with lfsr <
#'     `lfsr_thresh` in >= 2 / exactly 1 conditions), `classification`
#'     ("shared" if frac_top_shared >= 0.5, "private" if frac_top_private >=
#'     0.5, else "partial"). D-2a caveat: "shared" means active-in-both, not
#'     proven shared structure — at high loading SNR, independent loadings
#'     also reach frac_top_shared >= 0.5 (forced-pair diagnostic: 30/30 false
#'     "shared" at loading SE = 0.125 while mashr's mixture put 96% weight on
#'     the identity pattern, p_all_shared ~ 0.001). Check `p_all_shared`.
#'     Each per-group result also carries the full `mash` fit;
#'   * `summary`: one row per program. Programs in no matched group are
#'     classified "private" with evidence "no_significant_match".
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 20, n_genes = 150, program_size = 20, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 100)
#' fit <- canonicalize_programs(fit)
#' # Assess one matched set explicitly (skips the permutation matching step);
#' # with matches = NULL, candidate matches are data-driven.
#' if (requireNamespace("mashr", quietly = TRUE)) {
#'   grp <- data.frame(cell_type = c("CD8", "NK"), program_id = c("CD8_1", "NK_1"))
#'   res <- program_sharing(fit, matches = grp, n_top = 20)
#'   res$summary
#' }
#' @export
program_sharing <- function(fit, matches = NULL, alpha = 0.05, n_top = 50L,
                            lfsr_thresh = 0.05, n_perm = 199L, seed = 1L,
                            min_genes = 50L, min_null = 20L, z_thresh = 1,
                            data_driven = TRUE, npc = 3L, verbose = FALSE, ...) {
  .sharing_check_fit(fit, psd = TRUE)
  if (!.has_mashr()) {
    cli::cli_abort(c("Package {.pkg mashr} is required for {.fn program_sharing}.",
                     i = "Install it with {.code install.packages(\"mashr\")} and retry."))
  }
  match_table <- NULL
  if (is.null(matches)) {
    match_table <- .match_programs_loadings(fit, n_perm = n_perm, seed = seed)
    groups <- .sharing_groups(match_table, alpha = alpha)
  } else if (is.data.frame(matches) &&
             all(c("ct_a", "prog_a", "ct_b", "prog_b") %in% names(matches))) {
    match_table <- matches
    groups <- .sharing_groups(matches, alpha = alpha)
  } else if (is.data.frame(matches)) {
    groups <- list(.as_group(matches))
  } else if (is.list(matches)) {
    groups <- lapply(matches, .as_group)
  } else {
    .stopf("Unsupported `matches`; see ?program_sharing.")
  }
  .validate_groups(fit, groups)
  if (!is.null(matches)) {
    # D-2a: user-forced matches bypass the data-driven matching step, and the
    # "shared" label reflects the forced grouping. The classification counts
    # genes ACTIVE (lfsr < lfsr_thresh) in >= 2 conditions, which independent
    # loadings can also reach at high loading SNR (audit D-2a: 30/30 false
    # "shared" for forced independent pairs at loading SE = 0.125, while
    # mashr's own mixture put 96% weight on the identity pattern,
    # p_all_shared ~ 0.001). Check p_all_shared before interpreting "shared"
    # as shared structure.
    .warnf(paste0("`matches` supplied: matched sets are user-forced, not data-driven. ",
                   "The 'shared' classification counts top genes active (lfsr < %g) in >= 2 conditions ",
                   "and can label independent loadings 'shared' at high loading SNR; ",
                   "check `p_all_shared` (mixture weight of the equal-effects pattern) before interpretation."),
           lfsr_thresh)
  }

  results <- list()
  for (gi in seq_along(groups)) {
    g <- groups[[gi]]
    mats <- .sharing_matrices(fit, g, min_genes = min_genes)
    V <- .estimate_V(fit, g, z_thresh = z_thresh, min_null = min_null)
    data <- mashr::mash_set_data(mats$B, mats$S, V = V)
    U <- mashr::cov_canonical(data)
    if (data_driven) {
      U <- tryCatch({
        m1 <- mashr::mash_1by1(data)
        strong <- mashr::get_significant_results(m1, thresh = alpha)
        if (length(strong) >= max(10L, 2L * ncol(mats$B))) {
          Ud <- mashr::cov_pca(data, npc = min(npc, ncol(mats$B)), subset = strong)
          c(U, mashr::cov_ed(data, Ud, subset = strong))
        } else {
          U
        }
      }, error = function(e) {
        .warnf("Data-driven covariances failed (%s); using canonical only.", conditionMessage(e))
        U
      })
    }
    run <- function(Ulist) {
      if (verbose) mashr::mash(data, Ulist = Ulist, ...) else suppressMessages(mashr::mash(data, Ulist = Ulist, ...))
    }
    # Data-driven covariances can be numerically non-PSD when the number of
    # conditions is small (e.g. 4 blocks, npc = 3 leaves one ~0 eigenvalue that
    # rounds negative); fall back to canonical-only in that case.
    res <- tryCatch(run(U), error = function(e) {
      if (!data_driven) stop(e)
      .warnf("mash fit with data-driven covariances failed (%s); refitting with canonical covariances only.",
             conditionMessage(e))
      run(mashr::cov_canonical(data))
    })
    # Posterior accessors: CRAN mashr (<= 0.2.x) does not export get_pm/get_lfsr;
    # the $result fields are stable across versions, with accessor fallback.
    pm <- if ("get_pm" %in% getNamespaceExports("mashr"))
      get("get_pm", envir = asNamespace("mashr"))(res) else res$result$PosteriorMean
    lf <- if ("get_lfsr" %in% getNamespaceExports("mashr"))
      get("get_lfsr", envir = asNamespace("mashr"))(res) else res$result$lfsr
    pw <- res$posterior_weights
    ee <- grep("^equal_effects", colnames(pw))
    p_ee <- if (length(ee)) rowSums(pw[, ee, drop = FALSE]) else rep(0, nrow(pw))
    top <- order(-apply(abs(pm), 1L, max))[seq_len(min(n_top, nrow(pm)))]
    nsig <- rowSums(lf < lfsr_thresh)
    frac_shared <- mean(nsig[top] >= 2L)
    frac_private <- mean(nsig[top] == 1L)
    class_g <- if (frac_shared >= 0.5) "shared"
               else if (frac_private >= 0.5) "private"
               else "partial"
    pi <- res$fitted_g$pi
    pat <- tapply(pi, sub("\\.\\d+$", "", names(pi)), sum)
    results[[gi]] <- list(
      cell_types = g$cell_type, programs = g$program_id, genes = mats$genes,
      B_hat = mats$B, S_hat = mats$S, V = V,
      posterior_mean = pm, lfsr = lf,
      mixture_weights = pi, pattern_weights = pat,
      p_all_shared = mean(p_ee[top]),
      frac_top_shared = frac_shared, frac_top_private = frac_private,
      classification = class_g, mash = res
    )
  }
  # NB: paste0() treats zero-length args as "" (see canonicalize_programs), so
  # guard the no-groups case: paste0("group", integer(0)) would return "group".
  if (length(results)) names(results) <- paste0("group", seq_along(results))

  allp <- .program_table(fit)
  summ <- do.call(rbind, lapply(seq_len(nrow(allp)), function(i) {
    prog <- allp$program_id[i]; ct <- allp$cell_type[i]
    hit <- which(vapply(groups, function(g) any(g$program_id == prog & g$cell_type == ct),
                        logical(1L)))
    if (length(hit)) {
      r <- results[[hit[1L]]]
      data.frame(program_id = prog, cell_type = ct, group_id = hit[1L],
                 n_conditions = length(r$cell_types), classification = r$classification,
                 frac_top_shared = r$frac_top_shared,
                 frac_top_private = r$frac_top_private,
                 p_all_shared = r$p_all_shared, evidence = "mashr",
                 stringsAsFactors = FALSE)
    } else {
      data.frame(program_id = prog, cell_type = ct, group_id = NA_integer_,
                 n_conditions = 1L, classification = "private",
                 frac_top_shared = NA_real_, frac_top_private = NA_real_,
                 p_all_shared = NA_real_, evidence = "no_significant_match",
                 stringsAsFactors = FALSE)
    }
  }))
  rownames(summ) <- NULL
  out <- list(matches = match_table, groups = groups, results = results,
              summary = summ, alpha = alpha, n_top = n_top,
              lfsr_thresh = lfsr_thresh)
  class(out) <- "program_sharing"
  out
}

#' Print a program-sharing summary
#'
#' @param x A `program_sharing` object from [program_sharing()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.program_sharing <- function(x, ...) {
  cat(sprintf("program_sharing: %d matched group(s), %d program(s) assessed\n",
              length(x$results), nrow(x$summary)))
  for (i in seq_along(x$results)) {
    r <- x$results[[i]]
    cat(sprintf("  group%d [%s]: %s (frac shared %.2f / private %.2f, P(all-shared) %.2f)\n",
                i, paste(r$programs, collapse = " + "), r$classification,
                r$frac_top_shared, r$frac_top_private, r$p_all_shared))
  }
  n_un <- sum(is.na(x$summary$group_id))
  if (n_un) cat(sprintf("  %d unmatched program(s): private (no significant match)\n", n_un))
  invisible(x)
}
