## annot_theme.R -- Tier-2 theme annotation per CMP node (ULM, decoupler-style):
## for node loadings w (genes in the cell type's feature set), fit w ~ 1 + 1[g in set] for every set with
## >= min_size genes in the universe; t-statistic of the indicator; BH across sets within node.
## Returns per-node top themes. Gene-overlap class of a CMP (Tier-2): fraction of member-node pairs whose
## top-|loading| gene sets (top_n) overlap more than expected (hypergeometric, BH).
annotate_themes <- function(fit, cmp, sets, min_size = 10L, top = 3L, q = 0.05, use_abs = TRUE) {
  do.call(rbind, lapply(seq_len(nrow(cmp)), function(i) {
    ct <- cmp$cell_type[i]; pr <- cmp$program_id[i]
    w <- fit$loadings[[ct]][, pr]; g <- names(w); if (use_abs) w <- abs(w)
    res <- lapply(names(sets), function(s) {
      x <- as.numeric(g %in% sets[[s]]); if (sum(x) < min_size) return(NULL)
      m <- stats::lm.fit(cbind(1, x), w); rss <- sum(m$residuals^2); df <- length(w) - 2
      se <- sqrt(rss / df / sum((x - mean(x))^2)); t <- m$coefficients[2] / se
      data.frame(set = s, t = t, p = 2 * stats::pt(-abs(t), df))
    })
    r <- do.call(rbind, res); if (is.null(r)) return(NULL)
    r$padj <- stats::p.adjust(r$p, "BH"); r <- r[order(-abs(r$t)), ]
    r <- utils::head(r[r$padj < q, , drop = FALSE], top)
    if (!nrow(r)) return(data.frame(cluster = cmp$cluster[i], cell_type = ct, program_id = pr, set = NA, t = NA, padj = NA))
    data.frame(cluster = cmp$cluster[i], cell_type = ct, program_id = pr, set = r$set, t = r$t, padj = r$padj)
  }))
}
gene_overlap_class <- function(fit, cmp, top_n = 100L, q = 0.05) {
  if (nrow(cmp) < 2L) return("single")
  tops <- lapply(seq_len(nrow(cmp)), function(i) { w <- fit$loadings[[cmp$cell_type[i]]][, cmp$program_id[i]]; names(sort(-abs(w)))[seq_len(min(top_n, length(w)))] })
  univ <- lapply(seq_len(nrow(cmp)), function(i) rownames(fit$loadings[[cmp$cell_type[i]]]))
  ps <- c()
  for (i in seq_len(nrow(cmp) - 1L)) for (j in (i + 1L):nrow(cmp)) {
    U <- intersect(univ[[i]], univ[[j]]); a <- intersect(tops[[i]], U); b <- intersect(tops[[j]], U)
    k <- length(intersect(a, b)); ps <- c(ps, stats::phyper(k - 1, length(a), length(U) - length(a), length(b), lower.tail = FALSE))
  }
  f <- mean(stats::p.adjust(ps, "BH") < q)
  if (f >= 0.5) "shared-genes" else if (f > 0) "partial-overlap" else "distinct-genes"
}
## context association of a CMP: consensus score (mean of standardized member scores, sign-aligned to the
## first node) regressed on a donor-level variable; returns estimate, t, p.
cmp_context_association <- function(fit, cmp, context) {
  S <- lapply(seq_len(nrow(cmp)), function(i) { v <- fit$scores[[cmp$cell_type[i]]][, cmp$program_id[i]]; stats::setNames(as.numeric(scale(v)), names(v)) })
  ids <- sort(unique(unlist(lapply(S, names))))
  M <- sapply(S, function(v) v[ids]); if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  sg <- sign(stats::cor(M, M[, 1], use = "pairwise")); sg[!is.finite(sg)] <- 1
  cons <- rowMeans(sweep(M, 2, as.numeric(sg), `*`), na.rm = TRUE)
  x <- context[ids]; ok <- is.finite(cons) & !is.na(x)
  m <- summary(stats::lm(cons[ok] ~ x[ok]))$coefficients
  data.frame(estimate = m[2, 1], t = m[2, 3], p = m[2, 4], n = sum(ok))
}
