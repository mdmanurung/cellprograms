# Tier-2 annotations of coordinated multicellular programs. Each annotation has
# its own null and none of them affects whether a CMP is called.

.cmp_table <- function(cmp) if (inherits(cmp, "cmp_result")) cmp$cmps else cmp

.cmp_node_weight <- function(fit, ct, prog) {
  W <- fit$loadings[[ct]]
  if (is.null(W) || !prog %in% colnames(W)) return(NULL)
  w <- abs(W[, prog])
  lf <- fit$fits[[ct]]$flash$F_lfsr
  if (!is.null(lf) && nrow(lf) == length(w)) w <- w * (1 - lf[, match(prog, colnames(W))])
  w / max(w, 1e-12)
}

#' Ligand-receptor compatibility of coordinated programs
#'
#' For member programs a and b of a CMP in different cell types and a
#' ligand-receptor pair (L, R), scores `w_a[L] * w_b[R]`, where `w` is the
#' absolute gene loading times (1 - lfsr), scaled to a maximum of 1 within
#' each program. Significance uses a degree-preserving resource permutation
#' (ligand labels shuffled among ligands, receptor labels among receptors) with
#' the maximum over the CMP as the null statistic. Optionally, the donor-level
#' co-variation of the ligand in the sender and the receptor in the receiver is
#' reported from the pseudobulk (a MOFAtalk-style check).
#'
#' A hit means the coordinated programs are *compatible* with that interaction:
#' it is not evidence that the interaction causes the coordination. Programs
#' rich in secreted or surface genes (e.g. inflammatory programs) yield hits
#' more readily.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table.
#' @param lr Data frame with `ligand` and `receptor` columns (gene IDs matching
#'   the fit), e.g. a LIANA consensus resource.
#' @param n_perm Permutations (default 499).
#' @param alpha Level for reported pairs (default 0.05).
#' @param top Maximum pairs reported per CMP (default 10).
#' @param control `"fwer"` (default): each pair is compared with the permutation
#'   maximum over all pairs of the CMP (family-wise). `"fdr"`: each pair is
#'   compared with the pooled permutation null of all pair scores, then
#'   Benjamini-Hochberg over all tested (pair, sender, receiver) combinations;
#'   `p` is then the BH-adjusted value.
#' @param donor_check Add the donor-level ligand/receptor correlation (default TRUE).
#' @param seed Random seed.
#' @return Data frame: `cmp`, `sender`, `receiver`, `ligand`, `receptor`,
#'   `score`, `p` and, with `donor_check`, `donor_cor`, `donor_p` (Spearman
#'   test of ligand in sender vs receptor in receiver across shared donors) and
#'   `n_donors`.
#' @examples
#' # Toy canonicalized fit: program 1 of cell types A and B shares donor
#' # activity z (a planted coordinated program) and top genes G1-G15.
#' set.seed(1)
#' ids <- paste0("d", seq_len(40))
#' genes <- paste0("G", seq_len(600))
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
#' # plant a ligand at the top of A_1 and its receptor at the top of B_1
#' fit$loadings$A["G150", "A_1"] <- 5
#' fit$loadings$B["G151", "B_1"] <- 5
#' set.seed(2)
#' lr <- data.frame(ligand = c("G150", sample(paste0("G", 60:349), 300, replace = TRUE)),
#'                  receptor = c("G151", sample(paste0("G", 352:600), 300, replace = TRUE)))
#' annotate_cmp_lr(fit, res, lr, n_perm = 99)
#' @export
annotate_cmp_lr <- function(fit, cmp, lr, n_perm = 499L, alpha = 0.05, top = 10L,
                            control = c("fwer", "fdr"), donor_check = TRUE, seed = 1L) {
  control <- match.arg(control)
  tab <- .cmp_table(cmp)
  if (!all(c("ligand", "receptor") %in% names(lr))) .stopf("`lr` needs columns 'ligand' and 'receptor'.")
  lr <- unique(as.data.frame(lr)[, c("ligand", "receptor")])
  out <- lapply(unique(tab$cmp), function(k) {
    m <- tab[tab$cmp == k, , drop = FALSE]
    ws <- lapply(seq_len(nrow(m)), function(i) .cmp_node_weight(fit, m$cell_type[i], m$program_id[i]))
    names(ws) <- paste(m$cell_type, m$program_id, sep = "|")
    ws <- ws[!vapply(ws, is.null, TRUE)]
    if (length(ws) < 2L) return(NULL)
    nn <- names(ws); cts_n <- sub("\\|.*", "", nn)
    genes_all <- unique(unlist(lapply(ws, names)))
    Wm <- matrix(0, length(genes_all), length(ws), dimnames = list(genes_all, nn))
    for (j in seq_along(ws)) Wm[names(ws[[j]]), j] <- ws[[j]]
    ok <- outer(cts_n, cts_n, `!=`)                          # sender a, receiver b in different cell types
    gw <- function(g) { x <- Wm[match(g, genes_all), , drop = FALSE]; x[is.na(x)] <- 0; x }
    ## max over (a, b) of L[i, a] * R[i, b] for every pair i, computed one sender at a time
    max_score <- function(L, R) {
      m <- 0
      for (a in seq_along(nn)) {
        b <- which(ok[a, ]); if (!length(b)) next
        m <- max(m, L[, a] * R[, b, drop = FALSE])
      }
      m
    }
    withr::local_seed(seed + k)
    L0 <- gw(lr$ligand); R0 <- gw(lr$receptor)
    obs <- do.call(rbind, lapply(seq_along(nn), function(a) {
      b <- which(ok[a, ]); if (!length(b)) return(NULL)
      S <- L0[, a] * R0[, b, drop = FALSE]
      data.frame(i = rep(seq_len(nrow(lr)), length(b)), sender = nn[a], receiver = rep(nn[b], each = nrow(lr)),
                 score = as.vector(S), stringsAsFactors = FALSE)
    }))
    obs <- obs[obs$score > 0, , drop = FALSE]
    if (is.null(obs) || !nrow(obs)) return(NULL)
    ul <- unique(lr$ligand); ur <- unique(lr$receptor)
    Lu <- gw(ul); Ru <- gw(ur); li <- match(lr$ligand, ul); ri <- match(lr$receptor, ur)
    n_tests <- nrow(lr) * sum(ok)
    so <- sort(obs$score); cnt <- numeric(length(so))
    mx <- vapply(seq_len(n_perm), function(bb) {
      pl <- sample(length(ul)); pr <- sample(length(ur))   # label permutation: ligand ul[j] takes weights of ul[pl[j]]
      L <- Lu[pl[li], , drop = FALSE]; R <- Ru[pr[ri], , drop = FALSE]
      if (control == "fwer") return(max_score(L, R))
      for (a in seq_along(nn)) {                            # pooled null: count null scores >= each observed score
        b <- which(ok[a, ]); if (!length(b)) next
        v <- sort(as.vector(L[, a] * R[, b, drop = FALSE]))
        cnt <<- cnt + (length(v) - findInterval(so, v, left.open = TRUE))
      }
      0
    }, 1)
    if (control == "fwer") {
      obs$p <- vapply(obs$score, function(s) (1 + sum(mx >= s)) / (n_perm + 1), 1)
    } else {
      p <- (1 + cnt) / (1 + n_perm * n_tests)
      obs$p <- p[match(obs$score, so)]
      obs$p <- pmin(1, stats::p.adjust(c(obs$p, rep(1, n_tests - nrow(obs))), "BH")[seq_len(nrow(obs))])
    }
    obs <- obs[obs$p < alpha, , drop = FALSE]
    if (!nrow(obs)) return(NULL)
    obs$ligand <- lr$ligand[obs$i]; obs$receptor <- lr$receptor[obs$i]
    obs <- utils::head(obs[order(obs$p, -obs$score), ], top)
    data.frame(cmp = k, sender = obs$sender, receiver = obs$receiver, ligand = obs$ligand,
               receptor = obs$receptor, score = obs$score, p = obs$p, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    empty <- data.frame(cmp = integer(), sender = character(), receiver = character(), ligand = character(),
                        receptor = character(), score = numeric(), p = numeric(), stringsAsFactors = FALSE)
    if (donor_check) {
      # keep the documented schema stable even for an empty result
      empty$donor_cor <- numeric(); empty$donor_p <- numeric(); empty$n_donors <- numeric()
    }
    return(empty)
  }
  if (donor_check && is.null(fit$pseudobulk)) {
    .warnf("`fit` has no pseudobulk component; `donor_check` columns are NA.")
    out$donor_cor <- NA_real_; out$donor_p <- NA_real_; out$n_donors <- NA_real_
    donor_check <- FALSE
  }
  if (donor_check) {
    dc <- t(vapply(seq_len(nrow(out)), function(i) {
      cs <- sub("\\|.*", "", out$sender[i]); cr <- sub("\\|.*", "", out$receiver[i])
      Ys <- fit$pseudobulk[[cs]]; Yr <- fit$pseudobulk[[cr]]
      if (!out$ligand[i] %in% colnames(Ys) || !out$receptor[i] %in% colnames(Yr)) return(c(NA_real_, NA_real_, NA_real_))
      d <- intersect(rownames(Ys), rownames(Yr))
      if (length(d) < 5L) return(c(NA_real_, NA_real_, length(d)))
      ct <- suppressWarnings(stats::cor.test(Ys[d, out$ligand[i]], Yr[d, out$receptor[i]], method = "spearman", exact = FALSE))
      c(unname(ct$estimate), ct$p.value, length(d))
    }, numeric(3)))
    out$donor_cor <- dc[, 1]; out$donor_p <- dc[, 2]; out$n_donors <- dc[, 3]
  }
  rownames(out) <- NULL
  out
}

#' Themes of coordinated programs
#'
#' Univariate linear model (decoupler-style ULM) of each member program's
#' absolute loadings on gene-set membership, for every set with at least
#' `min_size` genes in the program's gene universe; BH across sets within a
#' program. Absolute loadings are used because program signs are arbitrary up
#' to canonicalization and a set can be enriched at either pole.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table.
#' @param sets Named list of gene vectors (e.g. MSigDB Hallmark).
#' @param min_size Minimum set size within the universe (default 10).
#' @param top Top sets per program (default 3).
#' @param q BH threshold (default 0.05).
#' @return Data frame: `cmp`, `cell_type`, `program_id`, `set`, `t`, `padj`.
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
#' sets <- list(planted = paste0("G", 1:15), other = paste0("G", 50:65))
#' annotate_cmp_themes(fit, res, sets, min_size = 5)
#' @export
annotate_cmp_themes <- function(fit, cmp, sets, min_size = 10L, top = 3L, q = 0.05) {
  tab <- .cmp_table(cmp)
  out <- lapply(seq_len(nrow(tab)), function(i) {
    ct <- tab$cell_type[i]; pr <- tab$program_id[i]
    w <- abs(fit$loadings[[ct]][, pr]); g <- names(w)
    r <- do.call(rbind, lapply(names(sets), function(s) {
      x <- as.numeric(g %in% sets[[s]])
      if (sum(x) < min_size || sum(x) == length(x)) return(NULL)
      m <- stats::lm.fit(cbind(1, x), w)
      df <- length(w) - 2
      se <- sqrt(sum(m$residuals^2) / df / sum((x - mean(x))^2))
      t <- unname(m$coefficients[2] / se)
      data.frame(set = s, t = t, p = 2 * stats::pt(-abs(t), df), stringsAsFactors = FALSE)
    }))
    none <- data.frame(cmp = tab$cmp[i], cell_type = ct, program_id = pr, set = NA_character_, t = NA_real_,
                       padj = NA_real_, stringsAsFactors = FALSE)
    if (is.null(r)) return(none)
    r$padj <- stats::p.adjust(r$p, "BH")
    r <- r[r$padj < q & r$t > 0, , drop = FALSE]
    r <- utils::head(r[order(-r$t), , drop = FALSE], top)
    if (!nrow(r)) return(none)
    data.frame(cmp = tab$cmp[i], cell_type = ct, program_id = pr, set = r$set, t = r$t, padj = r$padj,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Gene-overlap class of coordinated programs
#'
#' For each pair of member programs, tests whether their top-`top_n`
#' |loading| genes overlap more than expected within the shared gene universe
#' (hypergeometric, BH across pairs). A CMP is `"shared-genes"` if at least
#' half its pairs overlap, `"partial-overlap"` if some do and
#' `"distinct-genes"` otherwise. Distinct-gene CMPs are the case that
#' gene-matching methods cannot find.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table.
#' @param top_n Top genes per program (default 100).
#' @param q BH threshold (default 0.05).
#' @return Data frame: `cmp`, `n_pairs`, `frac_overlapping`, `class`.
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
#' cmp_gene_overlap(fit, res, top_n = 20)
#' @keywords internal
cmp_gene_overlap <- function(fit, cmp, top_n = 100L, q = 0.05) {
  tab <- .cmp_table(cmp)
  out <- lapply(unique(tab$cmp), function(k) {
    m <- tab[tab$cmp == k, , drop = FALSE]
    if (nrow(m) < 2L) {
      # no member pair to test; mean(logical(0)) would be NaN and crash the
      # class branch on an NA comparison
      return(data.frame(cmp = k, n_pairs = 0L, frac_overlapping = NA_real_,
                        class = "single-member", stringsAsFactors = FALSE))
    }
    W <- lapply(seq_len(nrow(m)), function(i) fit$loadings[[m$cell_type[i]]][, m$program_id[i]])
    tops <- lapply(W, function(w) names(sort(-abs(w)))[seq_len(min(top_n, length(w)))])
    ps <- c()
    for (i in seq_len(nrow(m) - 1L)) for (j in (i + 1L):nrow(m)) {
      U <- intersect(names(W[[i]]), names(W[[j]]))
      a <- intersect(tops[[i]], U); b <- intersect(tops[[j]], U)
      ps <- c(ps, stats::phyper(length(intersect(a, b)) - 1, length(a), length(U) - length(a), length(b),
                                lower.tail = FALSE))
    }
    f <- mean(stats::p.adjust(ps, "BH") < q)
    data.frame(cmp = k, n_pairs = length(ps), frac_overlapping = f,
               class = if (f >= 0.5) "shared-genes" else if (f > 0) "partial-overlap" else "distinct-genes",
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  if (is.null(res)) {
    # empty CMP table: return the documented schema, not NULL
    res <- data.frame(cmp = integer(0), n_pairs = integer(0),
                      frac_overlapping = numeric(0), class = character(0),
                      stringsAsFactors = FALSE)
  }
  res
}

#' Association of coordinated programs with a donor variable
#'
#' Regresses each CMP's consensus donor score (`cmp_scores()`) on a donor-level
#' variable, optionally adjusting for covariates.
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table.
#' @param variable A `sample_metadata` column name or a vector named by observation ID.
#' @param covariates Optional character vector of `sample_metadata` columns to adjust for.
#' @return Data frame: `cmp`, `term`, `estimate`, `t`, `p`, `n`.
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
#' cmp_context_association(fit, res, "group")
#' @keywords internal
cmp_context_association <- function(fit, cmp, variable, covariates = NULL) {
  S <- cmp_scores(fit, cmp)
  if (any(covariates %in% c("y", ".ctx"))) {
    cli::cli_abort("{.arg covariates} must not use the reserved names {.val y} or {.val .ctx}.")
  }
  v <- .cmp_resolve_var(fit, variable, "variable", allow_na = TRUE)   # donors with NA are dropped
  md <- fit$sample_metadata[match(rownames(S), fit$sample_metadata$observation_id), , drop = FALSE]
  out <- lapply(seq_len(ncol(S)), function(j) {
    d <- data.frame(y = S[, j], .ctx = v[rownames(S)])
    for (cv in covariates) d[[cv]] <- md[[cv]]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    cvs <- covariates[vapply(covariates, function(cv) length(unique(d[[cv]])) > 1L, TRUE)]   # drop constant covariates
    fm <- stats::as.formula(paste("y ~ .ctx", paste(c("", cvs), collapse = " + ")))
    co <- summary(stats::lm(fm, d))$coefficients
    # Exact context-term selection via the model-matrix term assignment (term 1
    # is .ctx). A prefix grep would also capture adjustment covariates whose
    # names start with the same letter (e.g. "xsite").
    mm <- stats::model.matrix(fm, d)
    ctx_terms <- colnames(mm)[attr(mm, "assign") == 1L]
    rr <- rownames(co) %in% ctx_terms
    data.frame(cmp = sub("^CMP", "", colnames(S)[j]), term = rownames(co)[rr], estimate = co[rr, 1],
               t = co[rr, 3], p = co[rr, 4], n = nrow(d), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
