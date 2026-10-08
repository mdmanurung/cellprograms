## integrate.R — optional stage-2 integration of stage-1 programs
##
## Methods:
##   none       — no integration (default)
##   ica        — multi-start ICA on the combined program-score matrix
##   pva        — principal-vector alignment (Level 4) shared axes
##   block      — Level 3 candidate: block-structured joint EBMF
##   ev_bidifac — Level 3 candidate: EV-BIDIFAC linked decomposition
##
## Integration never replaces stage-1 programs; it layers on top with full
## traceability (state -> program weights -> gene loadings).

#' Integrate cellular programs into higher-order states
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param method Integration method.
#' @param level For `block`/`ev_bidifac`: operate on program `"scores"`
#'   (default) or stacked `"expression"`.
#' @param n_states Number of higher-order states (ICA: components; PVA:
#'   maximum shared directions; block/ev_bidifac: rank caps).
#' @param n_runs ICA multi-start runs for stability estimation.
#' @param impute For ICA only: `"error"` (default; require complete
#'   observations) or `"mean"` (explicitly requested mean imputation).
#' @param reference PVA reference cell type (default: most programs).
#' @param seed Random seed.
#' @param ... Method-specific arguments.
#' @return The fit with `fit$integration` populated.
#' @export
integrate_programs <- function(fit,
                               method = c("none", "ica", "pva", "block", "ev_bidifac", "pmd"),
                               level = c("scores", "expression"),
                               n_states = NULL,
                               n_runs = 20L,
                               impute = c("error", "mean"),
                               reference = NULL,
                               seed = NULL,
                               ...) {
  method <- match.arg(method)
  level <- match.arg(level)
  if (!is.null(seed)) set.seed(seed)

  integ <- switch(
    method,
    none = list(method = "none"),
    ica  = .integrate_ica(fit, n_states, n_runs, match.arg(impute), seed),
    pva  = .integrate_pva(fit, n_states, reference, ...),
    block = fit_joint_block(fit, level = level, max_factors = n_states %||% 20L,
                            seed = seed, ...),
    ev_bidifac = fit_joint_evbidifac(fit, level = level, seed = seed, ...),
    pmd = .integrate_pmd(fit, n_states = n_states, seed = seed, ...)
  )
  fit$integration <- integ
  fit
}

## ---------------------------------------------------------------------------
## ICA on the combined program-score matrix
## ---------------------------------------------------------------------------

.integrate_ica <- function(fit, n_states, n_runs, impute, seed) {
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    .stopf("Package 'fastICA' required for method = 'ica'.")
  }
  Z <- program_scores(fit, "wide")
  cc <- stats::complete.cases(Z)
  if (!all(cc)) {
    if (impute == "error") {
      .stopf(paste0(
        "%d/%d observations have incomplete cell-type coverage. ",
        "ICA requires complete observations; rerun with impute = 'mean' ",
        "to explicitly request mean imputation."),
        sum(!cc), nrow(Z))
    }
    for (j in seq_len(ncol(Z))) {
      Z[is.na(Z[, j]), j] <- mean(Z[, j], na.rm = TRUE)
    }
    cc <- rep(TRUE, nrow(Z))
  }
  Zc <- Z[cc, , drop = FALSE]
  K <- n_states %||% min(10L, ncol(Zc) - 1L)
  K <- max(1L, min(K, ncol(Zc)))

  ## multi-start ICA with component matching (Icasso-style stability)
  runs <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    set.seed((seed %||% 0L) + r)
    runs[[r]] <- fastICA::fastICA(Zc, n.comp = K, alg.typ = "parallel",
                                  fun = "logcosh", verbose = FALSE,
                                  row.norm = FALSE, maxit = 500)
  }
  ## match components to run 1 by |correlation| of sample scores
  ref <- runs[[1L]]$S
  stab <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    cors <- vapply(runs[-1L], function(run) {
      max(abs(stats::cor(ref[, k], run$S)))
    }, numeric(1L))
    stab[k] <- mean(cors)
  }
  S <- runs[[1L]]$S; A <- runs[[1L]]$A
  rownames(S) <- rownames(Zc)
  state_ids <- sprintf("IC%02d", seq_len(K))
  colnames(S) <- state_ids
  ## A: programs x states mixing weights (fastICA: X = S A)
  prog_w <- t(A)
  rownames(prog_w) <- colnames(Zc); colnames(prog_w) <- state_ids

  list(
    method = "ica",
    states = S,
    program_weights = prog_w,
    stability = stats::setNames(stab, state_ids),
    n_runs = n_runs,
    complete_observations = rownames(Zc),
    imputed = impute == "mean" && !all(stats::complete.cases(program_scores(fit, "wide")))
  )
}

## ---------------------------------------------------------------------------
## PVA: principal-vector alignment shared axes
## ---------------------------------------------------------------------------

.integrate_pva <- function(fit, n_states, reference = NULL, n_perm = 200L,
                           dedup_cor = 0.9, ...) {
  ## all-pairs alignment (reference-anchored star alignment misses partial
  ## sharing that does not involve the reference cell type)
  cts <- names(fit$loadings)
  cts <- cts[vapply(fit$loadings, function(W) ncol(W) > 0L, logical(1L))]
  pairs <- utils::combn(cts, 2L, simplify = FALSE)
  al <- lapply(pairs, function(pr) principal_angles(fit, pr[1L], pr[2L],
                                                    space = "scores",
                                                    n_perm = n_perm, ...))
  Z <- program_scores(fit, "wide")

  ## collect all significant shared directions with their sample projections
  cands <- list()
  for (pa in al) {
    shared <- which(pa$p_values < 0.05 & pa$cosines > pa$asymptotic_ref)
    for (i in shared) {
      proj <- rep(NA_real_, nrow(Z)); names(proj) <- rownames(Z)
      n_contrib <- rep(0L, nrow(Z))
      for (side in c("a", "b")) {
        w <- pa$program_weights[[side]][, i]
        progs <- names(w)[names(w) %in% colnames(Z)]
        if (length(progs) == 0L) next
        pr <- as.numeric(Z[, progs, drop = FALSE] %*% w[progs])
        add <- !is.na(pr)
        cur <- proj[add]; cur[is.na(cur)] <- 0
        proj[add] <- (cur * n_contrib[add] + pr[add]) / (n_contrib[add] + 1)
        n_contrib[add] <- n_contrib[add] + 1L
      }
      cands[[length(cands) + 1L]] <- list(
        pair = c(pa$ct_a, pa$ct_b), direction = i,
        cosine = pa$cosines[i], p_value = pa$p_values[i],
        w_a = pa$program_weights$a[, i], w_b = pa$program_weights$b[, i],
        proj = proj
      )
    }
  }
  if (length(cands) == 0L) {
    return(list(method = "pva", states = NULL, axes = list(),
                note = "no shared directions above the permutation null"))
  }

  ## rank by cosine; greedily deduplicate axes with highly correlated
  ## sample projections (the same shared state seen through many pairs)
  ord <- order(-vapply(cands, `[[`, numeric(1L), "cosine"))
  cands <- cands[ord]
  kept <- list()
  for (cand in cands) {
    if (length(kept) > 0L) {
      cors <- vapply(kept, function(k) {
        abs(stats::cor(k$proj, cand$proj, use = "pairwise.complete.obs"))
      }, numeric(1L))
      if (any(cors > dedup_cor, na.rm = TRUE)) next
    }
    kept[[length(kept) + 1L]] <- cand
    if (!is.null(n_states) && length(kept) >= n_states) break
  }

  state_ids <- sprintf("PV%02d", seq_along(kept))
  S <- matrix(NA_real_, nrow(Z), length(kept),
              dimnames = list(rownames(Z), state_ids))
  for (j in seq_along(kept)) {
    S[, j] <- kept[[j]]$proj
    kept[[j]]$state_id <- state_ids[j]
    kept[[j]]$proj <- NULL
  }

  list(
    method = "pva",
    states = S,
    axes = kept,
    reference = reference,
    alignments = al
  )
}

## ---------------------------------------------------------------------------
## Traceability
## ---------------------------------------------------------------------------

#' Trace an integrated state back to genes
#'
#' For state with program weights a_j (programs j in cell type c), the
#' cell-type-specific gene summary is g_c = sum_j a_j w_j.
#'
#' @param fit A `cell_program_fit` with `fit$integration` set (or pass `integ`).
#' @param integ Optional integration result (defaults to `fit$integration`).
#' @param state State ID, e.g. `"IC03"` or `"PV01"`.
#' @param top_n Return only the top genes per cell type by |weight| (NULL = all).
#' @return A named list of per-cell-type gene weight vectors (sorted).
#' @export
trace_state <- function(fit, integ = NULL, state, top_n = NULL) {
  integ <- integ %||% fit$integration
  if (is.null(integ) || identical(integ$method, "none")) {
    .stopf("No integration present; run integrate_programs() first.")
  }
  pw <- switch(
    integ$method,
    ica = {
      if (!state %in% colnames(integ$program_weights)) .stopf("Unknown state '%s'.", state)
      integ$program_weights[, state]
    },
    pva = {
      ax <- Filter(function(a) identical(a$state_id, state), integ$axes)
      if (length(ax) == 0L) .stopf("Unknown state '%s'.", state)
      ## combine the pair's program weights
      a <- ax[[1L]]
      c(a$w_a, a$w_b)
    },
    block = ,
    ev_bidifac = {
      if (!is.null(integ$gene_loadings)) {
        ## expression-level: per-cell-type gene loadings of the state directly
        out <- lapply(integ$gene_loadings, function(W) {
          if (!state %in% colnames(W)) return(NULL)
          g <- W[, state]
          g <- g[abs(g) > 0]
          g <- g[order(-abs(g))]
          if (!is.null(top_n)) g <- utils::head(g, top_n)
          g
        })
        return(Filter(Negate(is.null), out))
      }
      if (!state %in% colnames(integ$program_weights)) .stopf("Unknown state '%s'.", state)
      integ$program_weights[, state]
    },
    .stopf("trace_state not implemented for method '%s'.", integ$method)
  )

  out <- list()
  for (ct in names(fit$loadings)) {
    W <- fit$loadings[[ct]]
    progs <- intersect(names(pw), colnames(W))
    if (length(progs) == 0L) next
    g <- as.numeric(W[, progs, drop = FALSE] %*% pw[progs])
    names(g) <- rownames(W)
    g <- g[abs(g) > 0]
    g <- g[order(-abs(g))]
    if (!is.null(top_n)) g <- utils::head(g, top_n)
    out[[ct]] <- g
  }
  out
}

#' Summarize an integration result
#' @param fit A `cell_program_fit` with integration.
#' @export
integration_summary <- function(fit) {
  integ <- fit$integration
  if (is.null(integ)) return(data.frame())
  switch(
    integ$method,
    none = data.frame(method = "none"),
    ica = data.frame(method = "ica", state = colnames(integ$states),
                     stability = integ$stability, stringsAsFactors = FALSE),
    pva = data.frame(
      method = "pva",
      state = vapply(integ$axes, `[[`, character(1L), "state_id"),
      pair = vapply(integ$axes, function(a) paste(a$pair, collapse = "-"), character(1L)),
      cosine = vapply(integ$axes, `[[`, numeric(1L), "cosine"),
      p_value = vapply(integ$axes, `[[`, numeric(1L), "p_value"),
      stringsAsFactors = FALSE
    ),
    block = ,
    ev_bidifac = integ$summary %||% data.frame(method = integ$method)
  )
}
