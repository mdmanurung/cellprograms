# Exact one-to-one maximum-weight assignment, with unmatched rows allowed.
# DP over column masks is small here (<=6 truth nodes; <=3 recurring clusters).
# Input row/column ID order resolves equal optima deterministically.
coord_assign <- function(w, threshold = .5) {
  nr <- nrow(w); nc <- ncol(w)
  out <- rep(NA_integer_, nr)
  if (!nr || !nc) return(out)
  if (nc > 15L) stop("Truth assignment supports at most 15 targets")
  states <- 2^nc
  value <- rep(-Inf, states); value[1] <- 0
  choices <- vector("list", states); choices[[1]] <- integer()
  for (i in seq_len(nr)) {
    nv <- rep(-Inf, states); ch <- vector("list", states)
    for (mask in which(is.finite(value)) - 1L) {
      options <- c(NA_integer_, which(is.finite(w[i, ]) & w[i, ] >= threshold))
      for (j in options) {
        if (is.na(j)) { nextmask <- mask; gain <- 0 }
        else {
          bit <- bitwShiftL(1L, j - 1L)
          if (bitwAnd(mask, bit)) next
          nextmask <- bitwOr(mask, bit); gain <- w[i, j]
        }
        v <- value[mask + 1L] + gain
        if (v > nv[nextmask + 1L]) {
          nv[nextmask + 1L] <- v
          ch[[nextmask + 1L]] <- c(choices[[mask + 1L]], j)
        }
      }
    }
    value <- nv; choices <- ch
  }
  choices[[which.max(value)]]
}

coord_truth_fit <- function(d, adjusted = TRUE) {
  scores <- d$truth$scores
  if (adjusted) for (ct in names(scores)) {
    scores[[ct]] <- .residualize(scores[[ct]], d$x$sample_metadata, c("institute", "nuisance_group"), ct)
    keep <- colSums(scores[[ct]]^2) > 1e-10
    scores[[ct]] <- scores[[ct]][, keep, drop = FALSE]
  }
  ld <- lapply(names(scores), function(ct) d$truth$loadings[[ct]][, colnames(scores[[ct]]), drop = FALSE])
  names(ld) <- names(scores)
  structure(list(scores = scores, loadings = ld, cell_types = names(scores),
    sample_metadata = d$x$sample_metadata,
    preprocessing = list(covariates = if (adjusted) c("institute", "nuisance_group") else NULL,
                         covariate_mode = "residualize")), class = "cell_program_fit")
}

coord_match_nodes <- function(fit, d, adjusted = TRUE) {
  reference <- coord_truth_fit(d, adjusted)
  rows <- list()
  for (ct in fit$cell_types) {
    Z <- fit$scores[[ct]]; W <- fit$loadings[[ct]]
    if (!ncol(Z)) next
    T <- reference$scores[[ct]][rownames(Z), , drop = FALSE]
    r <- suppressWarnings(abs(stats::cor(Z, T)))
    r[!is.finite(r)] <- 0
    r <- r[, order(colnames(r)), drop = FALSE]
    oi <- order(colnames(Z)); assigned <- coord_assign(r[oi, , drop = FALSE], .5)
    label <- rep(NA_character_, ncol(Z)); label[oi] <- colnames(r)[assigned]
    # Adjusted activity identity takes precedence over secondary gene identity.
    # Only an unmatched node can become the removed nuisance by loading cosine.
    if (adjusted && anyNA(label)) {
      g <- intersect(rownames(W), rownames(d$truth$loadings[[ct]]))
      nw <- d$truth$loadings[[ct]][g, "N"]
      norm <- sqrt(colSums(W[g, , drop = FALSE]^2)) * sqrt(sum(nw^2))
      co <- abs(as.vector(crossprod(W[g, , drop = FALSE], nw))) / pmax(norm, 1e-12)
      available <- which(is.na(label) & co >= .5)
      if (length(available)) label[available[order(-co[available], colnames(Z)[available])][1]] <- "N"
    }
    rows[[length(rows) + 1L]] <- data.frame(cell_type = ct, program_id = colnames(Z), truth = label)
  }
  if (length(rows)) do.call(rbind, rows) else data.frame(cell_type = character(), program_id = character(), truth = character())
}

coord_evaluate <- function(result, fit, d, method, adjusted = TRUE, oracle = FALSE) {
  memberships <- if (inherits(result, "coordinated_programs")) result$memberships else result$clusters
  breadth_only <- method == "pa_cc"
  map <- if (oracle) do.call(rbind, lapply(fit$cell_types, function(ct)
    data.frame(cell_type = ct, program_id = colnames(fit$scores[[ct]]), truth = colnames(fit$scores[[ct]])))) else coord_match_nodes(fit, d, adjusted)
  progs <- d$truth$programs
  targets <- progs[vapply(progs, function(p) isTRUE(p$activity_true), logical(1))]
  groups <- split(memberships, memberships$cluster)
  w <- matrix(0, length(groups), length(targets), dimnames = list(names(groups), names(targets)))
  labels <- lapply(groups, function(g) {
    idx <- match(paste(g$cell_type, g$program_id), paste(map$cell_type, map$program_id))
    map$truth[idx]
  })
  for (i in seq_along(groups)) for (j in seq_along(targets)) {
    g <- groups[[i]]; p <- targets[[j]]
    good <- g$cell_type %in% p$members
    if (!breadth_only) good <- good & !is.na(labels[[i]]) & labels[[i]] == p$id
    w[i, j] <- sum(good)
  }
  assigned <- coord_assign(w, 1)
  rows <- list(); correct <- 0L; found_clusters <- integer()
  for (nm in names(COORD_BREADTH)) {
    p <- targets[[nm]]
    if (is.null(p)) {
      rows[[nm]] <- data.frame(program = nm, present = FALSE, found = NA, recall = NA_real_, precision = NA_real_, breadth_error = NA_real_)
      next
    }
    j <- match(nm, names(targets)); i <- which(assigned == j)
    good <- if (length(i)) w[i, j] else 0L
    size <- if (length(i)) nrow(groups[[i]]) else 0L
    recall <- good / length(p$members); precision <- if (size) good / size else 0
    found <- recall >= .5 && precision >= .5 && good >= 2L
    if (found) found_clusters <- c(found_clusters, i)
    if (length(i)) correct <- correct + good
    rows[[nm]] <- data.frame(program = nm, present = TRUE, found = found, recall = recall, precision = precision,
                            breadth_error = size - length(p$members))
  }
  decoy <- if (breadth_only) NA else any(vapply(labels, function(v) sum(v == "N", na.rm = TRUE) >= 2L, logical(1)))
  # Clusters failing the fixed Found criterion, including duplicate fragments,
  # count as unmatched spurious outputs; nuisance is reported separately too.
  spurious <- any(!seq_along(groups) %in% found_clusters)
  fragments <- sum(vapply(seq_along(targets), function(j) sum(w[, j] >= 2L) > 1L, logical(1)))
  metric <- do.call(rbind, rows)
  metric$method <- method; metric$oracle <- oracle; metric$breadth_only <- breadth_only
  metric$spurious <- spurious; metric$decoy <- decoy; metric$n_clusters <- length(groups)
  metric$called_nodes <- if (breadth_only) NA_integer_ else nrow(memberships)
  metric$correct_nodes <- if (breadth_only) NA_integer_ else correct
  metric$fragments <- fragments
  metric
}
