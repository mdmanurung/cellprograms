## CMP stability under synchronized donor resampling (CP-2) -------------------
##
## assess_program_stability() subsamples observations independently per cell
## type, which is appropriate for local fit sensitivity but changes donor
## overlap arbitrarily between types; it cannot answer whether the *same
## donors* support an inferred coordinated multicellular program (CMP).
## assess_cmp_stability() draws donors ONCE per replicate, applies that draw
## to every cell type, and reruns CMP detection with the predeclared caller
## configuration. Subsampling is without replacement: duplicating donor rows
## would pseudo-replicate the donor-level evidence the permutation null
## evaluates (spec CP-2 gate: no pseudo-replication).
##
## Scope: this conditions on the fitted programs (score-level resampling) and
## reruns DETECTION; per-factor fit stability under resampling is answered by
## assess_program_stability(). The two compose: a stable CMP requires both.

#' Assess CMP stability under synchronized donor resampling
#'
#' Draws a donor subset once per replicate (optionally within donor-level
#' strata), applies the identical donor set to all cell types, and reruns
#' [find_coordinated_programs()] with the caller's configuration. Reports
#' per-CMP membership recovery (Jaccard), assignment ambiguity, sign
#' consistency, donor-score concordance, cell-type survival, and honest
#' failure denominators. Replicates that error are ledgered and excluded from
#' conditional rates; rates are `NA` when no valid replicate exists rather
#' than an artificial zero.
#'
#' @param fit A canonicalized fit (see [find_coordinated_programs()]).
#' @param donor_col Optional `sample_metadata` column giving the biological
#'   donor when it differs from `observation_ids`. All observations of a drawn
#'   donor are kept together in every cell type where present (cluster
#'   subsampling); missing cell types/visits are never synthesized.
#' @param visit_col Optional `sample_metadata` column identifying the
#'   visit/timepoint within donor. Required when `donor_col` is given and any
#'   donor contributes more than one observation: the estimand is then the
#'   donor-synchronized draw retaining all of a donor's visits, and the
#'   per-visit index is preserved in replicate metadata.
#' @param donor_strata Optional `sample_metadata` column with a donor-level
#'   stratification (must be constant within donor). The draw is taken within
#'   each stratum separately.
#' @param n_boot Number of resampling replicates (default 100).
#' @param subsample_frac Fraction of donors drawn per replicate, without
#'   replacement (default 0.8).
#' @param jaccard_threshold Minimum membership Jaccard between a reference CMP
#'   and its best-matching replicate CMP to count as recovered (default 0.5).
#' @param cmp_args Named list of arguments passed to
#'   [find_coordinated_programs()] unchanged (the predeclared caller
#'   configuration), e.g. `list(strata = "site", n_perm = 199, tau = 2,
#'   link_threshold = 3)`. `seed` is always overridden by the per-replicate
#'   seed recorded in the replicate log.
#' @param seed Base seed for reproducibility. Replicate `b` uses seed
#'   `seed + b - 1` for both the donor draw and the CMP permutation null.
#' @param progress Print progress every 10 replicates.
#'
#' @return An object of class `cmp_stability`, a list with
#'   \describe{
#'     \item{summary}{one row per reference CMP: `cmp`, `n_members`,
#'       `members`, `valid_refits`, `recovered_n`, `recovery_freq`
#'       (conditional on valid replicates; `NA` when none), `mean_jaccard`,
#'       `mean_ambiguity_margin` (best minus runner-up Jaccard),
#'       `sign_concordance` (fraction of recovered replicates whose member
#'       signs match the reference orientation up to one global flip),
#'       `donor_score_concordance` (mean correlation of CMP donor scores on
#'       shared donors among recovered replicates), `many_to_one_n`
#'       (replicates where this reference CMP shared its best-matching
#'       replicate CMP with another reference CMP).}
#'     \item{members}{one row per reference CMP member: recovery frequency in
#'       the matched replicate CMP and sign concordance given presence.}
#'     \item{cell_type_survival}{per cell type, fraction of valid replicates
#'       passing the `min_donors` filter.}
#'     \item{replicates}{per-replicate ledger: `replicate`, `seed`, `status`,
#'       `error_class`, `error_message`, `elapsed_seconds`, `n_donors_drawn`,
#'       `n_obs`, `n_cmps`, `eligible_rows` (per cell type), `surviving_cell_types`.}
#'     \item{reference}{the reference [find_coordinated_programs()] result.}
#'     \item{settings}{the arguments used.}
#'   }
#'   When the reference fit yields no CMPs, stability is undefined: the
#'   function returns an empty `summary`, zero-row `replicates`, and
#'   `status = "no_reference_cmps"` without spending replicates.
#'
#' @section Identity levels:
#' `observation_ids` are the score rows (observation level). `donor_col`
#' promotes resampling to the donor level; `visit_col` disambiguates repeated
#' observations per donor. Cohort/stratification effects on the null remain
#' the responsibility of `cmp_args$strata`.
#'
#' @examples
#' set.seed(1)
#' ids <- paste0("d", seq_len(40))
#' z <- rnorm(40)
#' mk <- function(ct) {
#'   S <- matrix(rnorm(40 * 2), 40, 2, dimnames = list(ids, paste0(ct, "_", 1:2)))
#'   S[, 1] <- z + rnorm(40, sd = 0.3)
#'   S
#' }
#' fit <- structure(list(
#'   scores = list(A = mk("A"), B = mk("B")),
#'   loadings = list(A = matrix(0, 5, 2), B = matrix(0, 5, 2)),
#'   cell_types = c("A", "B"), observation_ids = ids,
#'   sample_metadata = data.frame(observation_id = ids)),
#'   class = c("cell_program_fit", "cell_program_data", "list"))
#' stab <- assess_cmp_stability(fit, n_boot = 5, subsample_frac = 0.8,
#'                              cmp_args = list(n_perm = 99, ambient = "none",
#'                                              min_shared = 10, min_donors = 10))
#' stab$summary
#' @export
assess_cmp_stability <- function(fit, donor_col = NULL, visit_col = NULL,
                                 donor_strata = NULL, n_boot = 100L,
                                 subsample_frac = 0.8, jaccard_threshold = 0.5,
                                 cmp_args = list(), seed = NULL,
                                 progress = TRUE) {
  if (!is.list(cmp_args) || "fit" %in% names(cmp_args)) {
    .stopf("`cmp_args` must be a named list of find_coordinated_programs() arguments without `fit`.")
  }
  if (length(cmp_args) && (is.null(names(cmp_args)) || any(!nzchar(names(cmp_args))))) {
    .stopf("`cmp_args` must be fully named; unnamed elements would bind positionally in find_coordinated_programs().")
  }
  if (!is.numeric(subsample_frac) || length(subsample_frac) != 1L ||
      subsample_frac <= 0 || subsample_frac > 1) {
    .stopf("`subsample_frac` must be a single number in (0, 1].")
  }
  map <- .cmpstab_donor_map(fit, donor_col, visit_col, donor_strata)
  cmp_call <- c(list(fit = fit), cmp_args)
  cmp_call$seed <- if (is.null(seed)) 1L else seed
  ref <- do.call(find_coordinated_programs, cmp_call)
  settings <- list(donor_col = donor_col, visit_col = visit_col,
                   donor_strata = donor_strata, n_boot = n_boot,
                   subsample_frac = subsample_frac,
                   jaccard_threshold = jaccard_threshold,
                   cmp_args = cmp_args, seed = seed)
  empty <- list(summary = data.frame(), members = data.frame(),
                cell_type_survival = data.frame(), replicates = data.frame(),
                reference = ref, settings = settings, status = "no_reference_cmps")
  class(empty) <- c("cmp_stability", "list")
  if (is.null(ref$cmps) || !nrow(ref$cmps)) {
    message("reference fit yields no CMPs; stability undefined (no replicates run).")
    return(empty)
  }

  ref_members <- split(paste(ref$cmps$cell_type, ref$cmps$program_id, sep = "|"),
                       ref$cmps$cmp)
  ref_signs <- split(ref$cmps$sign, ref$cmps$cmp)
  ref_scores <- cmp_scores(fit, ref)
  min_donors <- if (is.null(cmp_args$min_donors)) 20L else cmp_args$min_donors

  # No master-seed call here: with a non-NULL seed every replicate reseeds
  # deterministically (seed + b - 1), so a function-level seed would be dead;
  # with seed = NULL the replicate seeds are drawn from the user's RNG stream.
  rep_seeds <- if (is.null(seed)) sample.int(.Machine$integer.max, n_boot) else seed + seq_len(n_boot) - 1L

  K <- length(ref_members)
  acc <- vector("list", K)           # per reference CMP replicate rows
  for (k in seq_len(K)) {
    acc[[k]] <- list(best = numeric(0), runner = numeric(0), matched = integer(0),
                     concordant = logical(0), score_cor = numeric(0),
                     many_to_one = logical(0))
  }
  memb_hit <- vector("list", K)      # per member presence in matched CMP
  memb_sign <- vector("list", K)     # per member sign concordance given presence
  for (k in seq_len(K)) {
    memb_hit[[k]] <- stats::setNames(vector("list", length(ref_members[[k]])), ref_members[[k]])
    memb_sign[[k]] <- stats::setNames(vector("list", length(ref_members[[k]])), ref_members[[k]])
  }
  ct_survive <- stats::setNames(vector("list", length(fit$scores)), names(fit$scores))
  repl_log <- vector("list", n_boot)
  valid_refits <- 0L
  failed_refits <- 0L

  for (b in seq_len(n_boot)) {
    if (progress && b %% 10L == 0L) message(sprintf("cmp stability replicate %d/%d", b, n_boot))
    t0 <- Sys.time()
    rs <- rep_seeds[b]
    rec <- tryCatch({
      withr::local_seed(rs)
      drawn <- .cmpstab_draw(map$donor, subsample_frac, map$strata)
      obs_keep <- map$observation_id[map$donor %in% drawn]
      fit_b <- .cmpstab_build_fit(fit, obs_keep)
      call_b <- c(list(fit = fit_b), cmp_args)
      call_b$seed <- rs
      res_b <- do.call(find_coordinated_programs, call_b)
      list(ok = TRUE, drawn = drawn, obs_keep = obs_keep, fit_b = fit_b, res_b = res_b)
    }, error = function(e) list(ok = FALSE, error = e))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!rec$ok) {
      failed_refits <- failed_refits + 1L
      repl_log[[b]] <- data.frame(replicate = b, seed = rs, status = "error",
                                  error_class = class(rec$error)[1],
                                  error_message = conditionMessage(rec$error),
                                  elapsed_seconds = elapsed, n_donors_drawn = NA_integer_,
                                  n_obs = NA_integer_, n_cmps = NA_integer_,
                                  eligible_rows = "", surviving_cell_types = "",
                                  stringsAsFactors = FALSE)
      next
    }
    valid_refits <- valid_refits + 1L
    res_b <- rec$res_b
    rep_cmps <- res_b$cmps
    rep_members <- if (is.null(rep_cmps) || !nrow(rep_cmps)) list() else
      split(paste(rep_cmps$cell_type, rep_cmps$program_id, sep = "|"), rep_cmps$cmp)
    rep_signs <- if (is.null(rep_cmps) || !nrow(rep_cmps)) list() else
      split(rep_cmps$sign, rep_cmps$cmp)
    rep_scores <- if (length(rep_members)) cmp_scores(rec$fit_b, res_b) else NULL

    elig <- vapply(rec$fit_b$scores, function(S) if (is.null(S)) 0L else nrow(S), 1L)
    surv <- names(rec$fit_b$scores)[elig >= min_donors]
    for (ct in names(ct_survive)) ct_survive[[ct]] <- c(ct_survive[[ct]], ct %in% surv)
    repl_log[[b]] <- data.frame(replicate = b, seed = rs, status = "ok",
                                error_class = "", error_message = "",
                                elapsed_seconds = elapsed,
                                n_donors_drawn = length(unique(rec$drawn)),
                                n_obs = length(rec$obs_keep),
                                n_cmps = length(rep_members),
                                eligible_rows = paste(sprintf("%s:%d", names(elig), elig), collapse = ";"),
                                surviving_cell_types = paste(surv, collapse = ";"),
                                stringsAsFactors = FALSE)

    # match replicate CMPs to each reference CMP by membership Jaccard
    matched_ids <- rep(NA_integer_, K)
    for (k in seq_len(K)) {
      mk <- ref_members[[k]]
      if (!length(rep_members)) {
        acc[[k]]$best <- c(acc[[k]]$best, 0); acc[[k]]$runner <- c(acc[[k]]$runner, 0)
        acc[[k]]$matched <- c(acc[[k]]$matched, NA_integer_)
        acc[[k]]$concordant <- c(acc[[k]]$concordant, NA)
        acc[[k]]$score_cor <- c(acc[[k]]$score_cor, NA_real_)
        for (m in mk) memb_hit[[k]][[m]] <- c(memb_hit[[k]][[m]], FALSE)
        next
      }
      jac <- vapply(rep_members, function(rm) .jaccard(mk, rm), 1)
      ord <- order(-jac)
      best <- jac[ord[1]]; runner <- if (length(ord) > 1L) jac[ord[2]] else 0
      j <- as.integer(names(jac)[ord[1]])
      matched_ids[k] <- j
      acc[[k]]$best <- c(acc[[k]]$best, best); acc[[k]]$runner <- c(acc[[k]]$runner, runner)
      acc[[k]]$matched <- c(acc[[k]]$matched, j)
      # sign concordance over shared members, up to one global flip
      shared <- intersect(mk, rep_members[[as.character(j)]])
      ps <- if (length(shared)) {
        ref_signs[[k]][match(shared, mk)] * rep_signs[[as.character(j)]][match(shared, rep_members[[as.character(j)]])]
      } else {
        numeric(0)
      }
      ps <- ps[!is.na(ps)]
      flip <- NA_real_
      if (length(ps)) {
        s <- sign(sum(ps))
        # a perfectly balanced sign sum has no majority orientation:
        # concordance is undefined (NA), not discordant
        if (s != 0) flip <- s
      }
      acc[[k]]$concordant <- c(acc[[k]]$concordant,
                               if (is.na(flip)) NA else all(ps == flip))
      # donor-score concordance on shared donors
      sc <- NA_real_
      if (best >= jaccard_threshold && !is.null(rep_scores) && paste0("CMP", j) %in% colnames(rep_scores)) {
        a <- ref_scores[rownames(ref_scores) %in% rec$obs_keep, k]
        bb <- rep_scores[, paste0("CMP", j)]
        common <- intersect(names(a)[!is.na(a)], names(bb)[!is.na(bb)])
        if (length(common) >= 4L) {
          va <- a[common]; vb <- bb[common]
          if (stats::sd(va) > 0 && stats::sd(vb) > 0) sc <- stats::cor(va, vb)
        }
      }
      acc[[k]]$score_cor <- c(acc[[k]]$score_cor, sc)
      # member-level presence and sign in the matched replicate CMP
      for (m in mk) {
        present <- m %in% rep_members[[as.character(j)]]
        memb_hit[[k]][[m]] <- c(memb_hit[[k]][[m]], present)
        if (present) {
          pr <- ref_signs[[k]][match(m, mk)] * rep_signs[[as.character(j)]][match(m, rep_members[[as.character(j)]])]
          memb_sign[[k]][[m]] <- c(memb_sign[[k]][[m]], is.finite(pr) && !is.na(flip) && pr == flip)
        }
      }
    }
    # many-to-one: a replicate CMP that is the best match for > 1 reference CMP
    tab <- table(matched_ids[!is.na(matched_ids)])
    clash <- as.integer(names(tab)[tab > 1L])
    for (k in seq_len(K)) {
      acc[[k]]$many_to_one <- c(acc[[k]]$many_to_one, !is.na(matched_ids[k]) && matched_ids[k] %in% clash)
    }
  }

  summary <- do.call(rbind, lapply(seq_len(K), function(k) {
    best <- acc[[k]]$best
    rec_flag <- best >= jaccard_threshold
    data.frame(cmp = k, n_members = length(ref_members[[k]]),
               members = paste(ref_members[[k]], collapse = ";"),
               valid_refits = valid_refits, failed_refits = failed_refits,
               recovered_n = sum(rec_flag),
               recovery_freq = if (valid_refits > 0L) sum(rec_flag) / valid_refits else NA_real_,
               mean_jaccard = if (length(best)) mean(best) else NA_real_,
               mean_ambiguity_margin = if (length(best)) mean(best - acc[[k]]$runner) else NA_real_,
               sign_concordance = if (any(rec_flag)) mean(acc[[k]]$concordant[rec_flag], na.rm = TRUE) else NA_real_,
               donor_score_concordance = if (any(rec_flag)) mean(acc[[k]]$score_cor[rec_flag], na.rm = TRUE) else NA_real_,
               many_to_one_n = sum(acc[[k]]$many_to_one),
               stringsAsFactors = FALSE)
  }))
  members <- do.call(rbind, lapply(seq_len(K), function(k) {
    mk <- ref_members[[k]]
    do.call(rbind, lapply(mk, function(m) {
      hits <- memb_hit[[k]][[m]]
      sg <- memb_sign[[k]][[m]]
      data.frame(cmp = k, member = m,
                 cell_type = sub("\\|.*", "", m), program_id = sub(".*\\|", "", m),
                 recovery_freq = if (length(hits)) mean(unlist(hits)) else NA_real_,
                 sign_concordance = if (length(sg)) mean(unlist(sg)) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
  cts <- do.call(rbind, lapply(names(ct_survive), function(ct) {
    v <- ct_survive[[ct]]
    data.frame(cell_type = ct,
               survival_freq = if (length(v)) mean(unlist(v)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  out <- list(summary = summary, members = members, cell_type_survival = cts,
              replicates = do.call(rbind, repl_log), reference = ref,
              settings = settings, status = "ok")
  class(out) <- c("cmp_stability", "list")
  out
}

#' @keywords internal
#' @export
print.cmp_stability <- function(x, ...) {
  cat("CMP stability (", x$settings$n_boot, " replicates, fraction ",
      x$settings$subsample_frac, ")\n", sep = "")
  if (identical(x$status, "no_reference_cmps")) {
    cat("reference fit yielded no CMPs; stability undefined.\n")
    return(invisible(x))
  }
  print(x$summary[, c("cmp", "n_members", "valid_refits", "recovered_n",
                      "recovery_freq", "mean_jaccard", "sign_concordance")],
        row.names = FALSE)
  invisible(x)
}

## ---- internals --------------------------------------------------------------

# Validated observation -> donor (+ visit, strata) map.
.cmpstab_donor_map <- function(fit, donor_col, visit_col, donor_strata) {
  ids <- fit$observation_ids
  if (is.null(ids)) .stopf("`fit$observation_ids` is required.")
  for (S in fit$scores) {
    if (is.null(S)) next
    if (is.null(rownames(S)) || !all(rownames(S) %in% ids)) {
      .stopf("score matrices must have rownames drawn from fit$observation_ids.")
    }
  }
  md <- fit$sample_metadata
  if (is.null(donor_col)) {
    return(data.frame(observation_id = ids, donor = ids, visit = NA_character_,
                      strata = NA_character_, stringsAsFactors = FALSE))
  }
  if (is.null(md) || !"observation_id" %in% names(md)) {
    .stopf("sample_metadata with an observation_id column is required when `donor_col` is used.")
  }
  if (!donor_col %in% names(md)) .stopf("donor_col '%s' not found in sample_metadata.", donor_col)
  miss <- setdiff(ids, md$observation_id)
  if (length(miss)) .stopf("sample_metadata is missing %d observation(s), e.g. '%s'.", length(miss), miss[1])
  md <- md[match(ids, md$observation_id), , drop = FALSE]
  donor <- as.character(md[[donor_col]])
  if (anyNA(donor)) .stopf("donor_col '%s' contains NA values.", donor_col)
  visit <- if (is.null(visit_col)) NA_character_ else {
    if (!visit_col %in% names(md)) .stopf("visit_col '%s' not found in sample_metadata.", visit_col)
    as.character(md[[visit_col]])
  }
  if (anyDuplicated(md$observation_id)) .stopf("sample_metadata observation_id must be unique.")
  n_per_donor <- table(donor)
  if (any(n_per_donor > 1L) && is.null(visit_col)) {
    .stopf(paste0("%d donors contribute > 1 observation; supply `visit_col` (or `observation_group_col`) ",
                  "so the paired-observation estimand is explicit."), sum(n_per_donor > 1L))
  }
  st <- if (is.null(donor_strata)) NA_character_ else {
    if (!donor_strata %in% names(md)) .stopf("donor_strata '%s' not found in sample_metadata.", donor_strata)
    sv <- as.character(md[[donor_strata]])
    bad <- tapply(sv, donor, function(x) length(unique(x)), simplify = TRUE)
    if (any(bad > 1L)) .stopf("donor_strata '%s' must be constant within donor.", donor_strata)
    sv
  }
  data.frame(observation_id = ids, donor = donor, visit = visit, strata = st,
             stringsAsFactors = FALSE)
}

# Draw donors once, without replacement, optionally within donor-level strata.
.cmpstab_draw <- function(donors, fraction, strata = NULL) {
  u <- unique(donors)
  if (is.null(strata) || all(is.na(strata))) {
    n <- max(1L, round(length(u) * fraction))
    return(sample(u, n))
  }
  st <- strata[match(u, donors)]
  unlist(lapply(split(u, st), function(g) {
    n <- max(1L, round(length(g) * fraction))
    sample(g, n)
  }), use.names = FALSE)
}

# Subset a canonicalized fit to the given observations (no synthesis, no
# duplication): score rows, observation_ids, and sample_metadata are filtered;
# loadings and program structure are unchanged (score-level resampling).
.cmpstab_build_fit <- function(fit, obs_keep) {
  obs_keep <- unique(obs_keep)
  ids_b <- fit$observation_ids[fit$observation_ids %in% obs_keep]
  fit$scores <- lapply(fit$scores, function(S) {
    if (is.null(S)) return(NULL)
    S[intersect(ids_b, rownames(S)), , drop = FALSE]
  })
  fit$observation_ids <- ids_b
  if (!is.null(fit$sample_metadata) && "observation_id" %in% names(fit$sample_metadata)) {
    fit$sample_metadata <- fit$sample_metadata[match(ids_b, fit$sample_metadata$observation_id), , drop = FALSE]
    rownames(fit$sample_metadata) <- NULL
  }
  fit
}

.jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) return(0)
  length(intersect(a, b)) / length(u)
}
