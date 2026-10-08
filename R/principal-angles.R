## principal-angles.R — Level 4: principal-angle program alignment
##
## Two geometries, answering different questions:
##
##   space = "scores"   : principal angles between the sample-activity
##                        subspaces Z_a and Z_b (over shared observations).
##                        Detects ACTIVITY COORDINATION — programs that covary
##                        across individuals regardless of which genes they
##                        use. This is the primary cross-cell-type sharing
##                        diagnostic (the multicellular-program notion).
##                        Null: permute sample labels of one cell type.
##
##   space = "loadings" : principal angles between the gene-loading subspaces
##                        W_a and W_b (over the union gene set). Detects
##                        SHARED GENE USAGE. Used for program matching across
##                        bootstrap refits / cohorts, and flags Scenario-7-type
##                        structure (same genes, possibly independent activity).
##                        Null: permute gene labels of one cell type.
##
## Math: orthonormalize to U_a, U_b; SVD of the cross-Gram U_a' U_b = P S Q'
## gives cosines sigma_i = cos(theta_i) and principal vectors U_a P, U_b Q.
## Rotation-invariant within each subspace by construction.
##
## Gene-universe policy (loadings space): intersection. Zero-padding to the
## union preserves CROSS inner products but NOT self-norms (each side's
## private genes inflate its own normalization), so padded angles are diluted
## by private genes — they are not equivalent to intersect-only angles.
## Intersection-only is the standard "shared gene usage" geometry; the
## intersection size is reported as the effective dimension for null
## calibration. (Identical gene panels => intersection == union, no effect.)

## ---------------------------------------------------------------------------
## generic engine: rows of Ma/Mb are the common space (samples or genes)
## ---------------------------------------------------------------------------

.pad_rows <- function(M, ids) {
  out <- matrix(0, nrow = length(ids), ncol = ncol(M),
                dimnames = list(ids, colnames(M)))
  common <- intersect(rownames(M), ids)
  out[common, ] <- M[common, , drop = FALSE]
  out
}

.subspace <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  r <- sum(s$d > tol * max(s$d))
  if (r == 0L) {
    return(list(U = matrix(0, nrow(M), 0L), Vsinvi = matrix(0, ncol(M), 0L), rank = 0L))
  }
  U <- s$u[, seq_len(r), drop = FALSE]
  Vsinvi <- s$v[, seq_len(r), drop = FALSE] %*% diag(1 / s$d[seq_len(r)], nrow = r)
  list(U = U, Vsinvi = Vsinvi, rank = r)
}

.pa_engine <- function(Ma, Mb, n_perm, seed) {
  sa <- .subspace(Ma); sb <- .subspace(Mb)
  if (sa$rank == 0L || sb$rank == 0L) .stopf("Degenerate (rank-0) subspace.")
  s <- svd(crossprod(sa$U, sb$U))
  sigma <- pmin(pmax(s$d, 0), 1)

  null_sigma <- NULL; p_values <- NULL; p_pair <- NA_real_
  z_pair <- excess_frac <- NA_real_
  if (n_perm > 0L) {
    if (!is.null(seed)) set.seed(seed)
    null_sigma <- matrix(NA_real_, nrow = n_perm, ncol = length(sigma))
    for (p in seq_len(n_perm)) {
      Mb_p <- Mb[sample.int(nrow(Mb)), , drop = FALSE]
      sb_p <- .subspace(Mb_p)
      d_p <- svd(crossprod(sa$U, sb_p$U))$d
      null_sigma[p, seq_along(d_p)] <- pmin(d_p, 1)
    }
    p_values <- vapply(seq_along(sigma), function(i) {
      (1 + sum(null_sigma[, i] >= sigma[i], na.rm = TRUE)) / (n_perm + 1)
    }, numeric(1L))
    ## one pair-level test: T = sum(cos^2) = ||Ua'Ub||_F^2 (rotation-invariant;
    ## no multiplicity over directions), same permutation null
    t_null <- rowSums(null_sigma^2, na.rm = TRUE)
    p_pair <- (1 + sum(t_null >= sum(sigma^2))) / (n_perm + 1)
    ## graded effect size, informative where p saturates at 1/(n_perm+1):
    ## z = (T - mean null) / sd null; excess_frac = (T - E0) / (min(ka,kb) - E0)
    ## is 0 at chance and 1 when the smaller subspace is fully contained
    t_obs <- sum(sigma^2); e0 <- mean(t_null)
    sd0 <- stats::sd(t_null)
    z_pair <- if (is.finite(sd0) && sd0 > 0) (t_obs - e0) / sd0 else NA_real_   # degenerate null: no z
    excess_frac <- (t_obs - e0) / (min(sa$rank, sb$rank) - e0)
  }

  list(
    cosines = sigma,
    angles = acos(sigma),
    p_values = p_values,
    p_pair = p_pair, z_pair = z_pair, excess_frac = excess_frac,
    null_mean = if (!is.null(null_sigma)) colMeans(null_sigma, na.rm = TRUE) else NULL,
    null_q95 = if (!is.null(null_sigma)) apply(null_sigma, 2L, stats::quantile, 0.95, na.rm = TRUE) else NULL,
    n_shared_05 = if (!is.null(p_values)) sum(p_values < 0.05) else NA_integer_,
    principal_vectors = list(a = sa$U %*% s$u, b = sb$U %*% s$v),
    program_weights = list(a = sa$Vsinvi %*% s$u, b = sb$Vsinvi %*% s$v),
    ranks = c(a = sa$rank, b = sb$rank)
  )
}

## residual-SD gene weighting (loadings space only). Adapter: the repo keeps
## per-gene residuals implicitly as Y_used - Z W' (canonical scores/loadings
## reproduce the reconstruction exactly), so derive the SD from that.
.gene_weights <- function(fit, ct_a, ct_b, genes) {
  rsd <- function(ct) {
    Y <- fit$fits[[ct]]$Y_used
    apply(Y - fit$scores[[ct]] %*% t(fit$loadings[[ct]]), 2L, stats::sd)
  }
  sd_a <- rsd(ct_a); sd_b <- rsd(ct_b)
  wa <- rep(NA_real_, length(genes)); names(wa) <- genes
  wb <- wa
  wa[intersect(genes, names(sd_a))] <- sd_a[intersect(genes, names(sd_a))]
  wb[intersect(genes, names(sd_b))] <- sd_b[intersect(genes, names(sd_b))]
  s <- pmin(ifelse(is.na(wa), Inf, wa), ifelse(is.na(wb), Inf, wb))
  s[!is.finite(s)] <- stats::median(s[is.finite(s)], na.rm = TRUE)
  s[s < 1e-8] <- 1e-8
  w <- 1 / s
  pmin(w, 10 * stats::median(w))   # cap: one near-perfectly fitted gene must not dominate the angles
}

#' Principal angles between the program subspaces of two cell types
#'
#' @param fit A `cell_program_fit` object.
#' @param ct_a,ct_b Cell-type labels.
#' @param space `"scores"` (activity coordination; default) or `"loadings"`
#'   (shared gene usage).
#' @param weight Loadings space only: `"none"` or `"residual_sd"` whitening.
#' @param n_perm Permutations for the label-permutation null (0 skips).
#' @param seed Random seed for permutations.
#' @return An object of class `principal_angles`.
#' @export
principal_angles <- function(fit, ct_a, ct_b,
                             space = c("scores", "loadings"),
                             weight = c("none", "residual_sd"),
                             n_perm = 200L, seed = NULL) {
  space <- match.arg(space)
  weight <- match.arg(weight)

  if (space == "scores") {
    Za <- fit$scores[[ct_a]]; Zb <- fit$scores[[ct_b]]
    if (is.null(Za) || is.null(Zb)) .stopf("Unknown cell type(s).")
    if (ncol(Za) == 0L || ncol(Zb) == 0L) {
      .stopf("Both cell types must retain >= 1 program.")
    }
    ids <- intersect(rownames(Za), rownames(Zb))
    if (length(ids) < 20L) {
      .warnf("Only %d shared observations for %s vs %s; angles may be unstable.",
             length(ids), ct_a, ct_b)
    }
    Ma <- Za[ids, , drop = FALSE]; Mb <- Zb[ids, , drop = FALSE]
    eff_dim <- length(ids)
    if (weight != "none") {
      .warnf("`weight` applies to loadings space only; ignored for scores.")
    }
  } else {
    W_a <- fit$loadings[[ct_a]]; W_b <- fit$loadings[[ct_b]]
    if (is.null(W_a) || is.null(W_b)) .stopf("Unknown cell type(s).")
    if (ncol(W_a) == 0L || ncol(W_b) == 0L) {
      .stopf("Both cell types must retain >= 1 program.")
    }
    genes <- intersect(rownames(W_a), rownames(W_b))
    g_inter <- length(genes)
    if (g_inter < 200L) {
      .warnf("Gene intersection for %s vs %s is %d (< 200); angles may be unstable.",
             ct_a, ct_b, g_inter)
    }
    if (g_inter < 10L) .stopf("Gene intersection for %s vs %s is %d (< 10).",
                              ct_a, ct_b, g_inter)
    Ma <- W_a[genes, , drop = FALSE]; Mb <- W_b[genes, , drop = FALSE]
    if (weight == "residual_sd") {
      gw <- .gene_weights(fit, ct_a, ct_b, genes)
      Ma <- Ma * gw; Mb <- Mb * gw
    }
    eff_dim <- g_inter
  }

  eng <- .pa_engine(Ma, Mb, n_perm, seed)
  k_a <- eng$ranks["a"]; k_b <- eng$ranks["b"]
  asymptotic_ref <- (sqrt(k_a) + sqrt(k_b)) / sqrt(max(eff_dim, 1))
  ## sharing requires BOTH beating the permutation null AND a minimum effect
  ## size (the random-subspace chance level); p-only flags near-orthogonal
  ## directions whenever the empirical null happens to be tight
  if (!is.null(eng$p_values)) {
    eng$n_shared_05 <- sum(eng$p_values < 0.05 & eng$cosines > asymptotic_ref)
  }

  rownames(eng$program_weights$a) <- colnames(Ma)
  rownames(eng$program_weights$b) <- colnames(Mb)

  out <- list(
    ct_a = ct_a, ct_b = ct_b, space = space,
    cosines = eng$cosines, angles = eng$angles,
    p_values = eng$p_values,
    p_pair = eng$p_pair, z_pair = eng$z_pair, excess_frac = eng$excess_frac,
    underpowered = unname(asymptotic_ref >= 1),
    null_mean = eng$null_mean, null_q95 = eng$null_q95,
    asymptotic_ref = asymptotic_ref,
    n_shared_05 = eng$n_shared_05,
    principal_vectors = eng$principal_vectors,
    program_weights = eng$program_weights,
    dims = c(k_a = unname(k_a), k_b = unname(k_b), effective_dim = eff_dim),
    weight = if (space == "loadings") weight else "none",
    n_perm = n_perm
  )
  class(out) <- "principal_angles"
  out
}

#' @export
print.principal_angles <- function(x, ...) {
  cat(sprintf("principal_angles (%s space): %s vs %s\n", x$space, x$ct_a, x$ct_b))
  cat(sprintf("  ranks %d x %d, effective dim %d\n",
              x$dims["k_a"], x$dims["k_b"], x$dims["effective_dim"]))
  cat("  cosines:", paste(round(x$cosines, 3), collapse = ", "), "\n")
  if (!is.null(x$p_values)) {
    cat("  perm p :", paste(round(x$p_values, 3), collapse = ", "), "\n")
    cat(sprintf("  shared directions (p < 0.05): %d  [asymptotic null ref %.3f]\n",
                x$n_shared_05, x$asymptotic_ref))
    cat(sprintf("  effect: z = %.1f, excess fraction = %.2f (0 chance, 1 nested)\n", x$z_pair, x$excess_frac))
    cat(sprintf("  pair-level p (sum cos^2): %.4f%s\n", x$p_pair,
                if (isTRUE(x$underpowered)) "  [UNDERPOWERED: chance cosine >= 1, not callable]" else ""))
  }
  invisible(x)
}

# Deterministic per-pair seed from a master seed and the pair's cell-type names.
.pair_seed <- function(seed, pair) {
  h <- sum(utf8ToInt(paste(pair, collapse = "\r")) * seq_len(nchar(paste(pair, collapse = "\r"))))
  as.integer((as.numeric(seed) + h) %% .Machine$integer.max)
}

#' Pairwise sharing spectra across all cell types
#'
#' @param fit A `cell_program_fit` object.
#' @param pairs Optional list of length-2 cell-type vectors; default all pairs.
#' @param space Passed to `principal_angles` (default `"scores"`).
#' @param n_perm Permutations per pair. Default (NULL) scales with the number of
#'   pairs, `max(999, 40 * #pairs)`, so that the permutation floor `1/(n_perm+1)`
#'   is below the BH level of a single true pair (`0.05 / #pairs`); a smaller
#'   value warns. Pairs get distinct seeds derived from `seed`.
#' @param ... Passed to `principal_angles`.
#' @return Per-pair `principal_angles` objects plus a summary table. `p_pair` is
#'   a single permutation p-value per pair (statistic sum of squared cosines);
#'   `z_pair` and `excess_frac` are graded effect sizes against the
#'   permutation null (z-score of the statistic; fraction of the smaller
#'   subspace's dimensions shared beyond chance); `q_pair` is BH across callable pairs; `underpowered` marks pairs whose
#'   chance-level cosine `(sqrt(ka)+sqrt(kb))/sqrt(n)` is >= 1 (excluded from
#'   BH, never called); `shared` = `q_pair < 0.05`. Needs `n_perm > 0`.
#' @export
sharing_spectrum <- function(fit, pairs = NULL, space = c("scores", "loadings"),
                             n_perm = NULL, ...) {
  space <- match.arg(space)
  cts <- names(fit$loadings)
  cts <- cts[vapply(fit$loadings, function(W) ncol(W) > 0L, logical(1L))]
  if (is.null(pairs)) pairs <- utils::combn(cts, 2L, simplify = FALSE)
  if (is.null(n_perm)) n_perm <- max(999L, 40L * length(pairs))
  if (n_perm > 0L && 1 / (n_perm + 1) > 0.05 / length(pairs)) {
    .warnf("n_perm = %d gives a permutation floor of %.4f, above the BH level of one true pair among %d (%.4f); BH can only call pairs if many hit the floor. Use n_perm >= %d.",
           n_perm, 1 / (n_perm + 1), length(pairs), 0.05 / length(pairs), 20L * length(pairs))
  }
  dots <- list(...)
  res <- lapply(pairs, function(pr) {
    a <- dots
    if (!is.null(a$seed)) a$seed <- .pair_seed(a$seed, pr)   # distinct stream per pair
    do.call(principal_angles, c(list(fit, pr[1L], pr[2L], space = space, n_perm = n_perm), a))
  })
  names(res) <- vapply(pairs, paste, character(1L), collapse = " vs ")

  summ <- do.call(rbind, lapply(res, function(r) {
    data.frame(pair = paste(r$ct_a, "vs", r$ct_b),
               k_a = unname(r$dims["k_a"]), k_b = unname(r$dims["k_b"]),
               effective_dim = unname(r$dims["effective_dim"]),
               min_angle_deg = round(min(r$angles) * 180 / pi, 1),
               n_shared_05 = r$n_shared_05,
               p_pair = r$p_pair, z_pair = r$z_pair, excess_frac = r$excess_frac,
               underpowered = r$underpowered,
               stringsAsFactors = FALSE)
  }))
  ## BH over the pairs that can be called; underpowered pairs (chance-level
  ## cosine >= 1) get q = NA and are never called shared
  summ$q_pair <- NA_real_
  ok <- !summ$underpowered & !is.na(summ$p_pair)
  summ$q_pair[ok] <- stats::p.adjust(summ$p_pair[ok], "BH")
  summ$shared <- ok & summ$q_pair < 0.05 & !is.na(summ$q_pair)
  rownames(summ) <- NULL
  list(pairs = res, summary = summ, space = space)
}

#' Align all cell types' program subspaces to a reference cell type
#'
#' @param fit A `cell_program_fit` object.
#' @param reference Reference cell type (default: most programs).
#' @param space Passed to `principal_angles` (default `"scores"`).
#' @param ... Passed to `principal_angles`.
#' @export
align_programs <- function(fit, reference = NULL, space = c("scores", "loadings"), ...) {
  space <- match.arg(space)
  cts <- names(fit$loadings)
  cts <- cts[vapply(fit$loadings, function(W) ncol(W) > 0L, logical(1L))]
  if (is.null(reference)) {
    reference <- cts[which.max(vapply(fit$loadings[cts], ncol, integer(1L)))]
  }
  others <- setdiff(cts, reference)
  out <- lapply(others, function(ct) principal_angles(fit, reference, ct, space = space, ...))
  names(out) <- others
  attr(out, "reference") <- reference
  attr(out, "space") <- space
  out
}

#' Match programs across two fits (bootstrap refits or cohorts)
#'
#' @param fit_a,fit_b `cell_program_fit` objects.
#' @param cell_type Cell type to match (must exist in both fits).
#' @param method `"subspace"` (default) or `"correlation"`. Factor-level
#'   `matches` are the same greedy loading-cosine matches in both; `"subspace"`
#'   additionally returns principal `angles` and a `split_map` (heuristic
#'   thresholds: cosine > 0.5, |weight| > 0.3; weights carry a singular-value
#'   scale and are not comparable across cell types).
#' @param space Subspace geometry for matching: `"loadings"` (default; gene
#'   identity) or `"scores"` (activity; requires shared observations).
#' @param sim_threshold Minimum |cosine| to report a factor-level match.
#' @export
match_programs <- function(fit_a, fit_b, cell_type,
                           method = c("subspace", "correlation"),
                           space = c("loadings", "scores"),
                           sim_threshold = 0.5) {
  method <- match.arg(method)
  space <- match.arg(space)
  W_a <- fit_a$loadings[[cell_type]]; W_b <- fit_b$loadings[[cell_type]]
  if (is.null(W_a) || is.null(W_b)) .stopf("Cell type '%s' missing from a fit.", cell_type)
  if (ncol(W_a) == 0L || ncol(W_b) == 0L) {
    return(list(matches = data.frame(), angles = NULL,
                note = "no programs on at least one side"))
  }

  genes <- intersect(rownames(W_a), rownames(W_b))
  if (length(genes) < 10L) {
    .stopf("Cell type '%s': gene intersection between fits is %d (< 10).",
           cell_type, length(genes))
  }
  Wa <- W_a[genes, , drop = FALSE]; Wb <- W_b[genes, , drop = FALSE]

  na <- sqrt(colSums(Wa^2)); nb <- sqrt(colSums(Wb^2))
  na[na == 0] <- 1; nb[nb == 0] <- 1
  cos_sim <- crossprod(Wa, Wb) / outer(na, nb)
  gm <- .greedy_match(cos_sim)
  matches <- data.frame(
    program_a = colnames(Wa)[gm$a],
    program_b = colnames(Wb)[gm$b],
    cosine = vapply(seq_len(nrow(gm)), function(i) cos_sim[gm$a[i], gm$b[i]], numeric(1L)),
    stringsAsFactors = FALSE
  )
  matches <- matches[abs(matches$cosine) >= sim_threshold, , drop = FALSE]

  angles <- NULL; split_map <- NULL
  if (method == "subspace") {
    Ma <- if (space == "loadings") Wa else {
      ids <- intersect(rownames(fit_a$scores[[cell_type]]), rownames(fit_b$scores[[cell_type]]))
      fit_a$scores[[cell_type]][ids, , drop = FALSE]
    }
    Mb <- if (space == "loadings") Wb else {
      ids <- intersect(rownames(fit_a$scores[[cell_type]]), rownames(fit_b$scores[[cell_type]]))
      fit_b$scores[[cell_type]][ids, , drop = FALSE]
    }
    sa <- .subspace(Ma); sb <- .subspace(Mb)
    s <- svd(crossprod(sa$U, sb$U))
    sigma <- pmin(s$d, 1)
    angles <- acos(sigma)
    A_a <- sa$Vsinvi %*% s$u; A_b <- sb$Vsinvi %*% s$v
    rownames(A_a) <- colnames(Ma); rownames(A_b) <- colnames(Mb)
    shared_dirs <- which(sigma > 0.5)
    split_map <- lapply(shared_dirs, function(i) {
      list(direction = i, cosine = sigma[i],
           programs_a = rownames(A_a)[abs(A_a[, i]) > 0.3],
           programs_b = rownames(A_b)[abs(A_b[, i]) > 0.3])
    })
  }

  list(matches = matches, angles = angles, split_map = split_map,
       method = method, space = space, cell_type = cell_type)
}
