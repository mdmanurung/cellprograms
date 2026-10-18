# cellprograms: visualization suite
#
# ggplot2-based plot functions, one per analysis stage in the tutorials.
# Every function returns a ggplot object (printable, ggsave-able) and
# degrades gracefully on degenerate inputs (no programs, no edges): it
# returns an annotated empty plot rather than erroring, so tutorial and
# report pipelines never break on a sparse fit.

## ---------------------------------------------------------------------------
## internal helpers
## ---------------------------------------------------------------------------

.theme_programs <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

## Empty-plot guard: tutorials call these on real fits where a cell type may
## have kept zero programs or a CMP search may have found zero edges.
.empty_plot <- function(msg) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 4) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11))
  message(msg)
  p
}

## Symmetric diverging palette for score-like quantities (colorblind-safe).
.diverging_pal <- function() c("#2166AC", "#F7F7F7", "#B2182B")

## Long-format melt of a dimnamed matrix: columns rows / cols / value.
.melt_matrix <- function(M) {
  d <- as.data.frame(as.table(as.matrix(M)), stringsAsFactors = FALSE)
  names(d) <- c("rows", "cols", "value")
  d
}

## Clip to symmetric limits (avoids a scales dependency for oob handling).
.clip_sym <- function(x, lim) pmax(pmin(x, lim), -lim)

## ---------------------------------------------------------------------------
## 1. program scores heatmap
## ---------------------------------------------------------------------------

#' Plot program scores across donors
#'
#' Tile heatmap of the canonicalized program score matrix (donors x programs),
#' optionally with donor-annotation strips (e.g. severity, batch) aligned
#' underneath.
#'
#' @param fit A canonicalized `cell_program_fit` object.
#' @param cell_types Cell types to include (default: all with >= 1 program).
#' @param annotation Optional character vector of column names in
#'   `fit$sample_metadata` to display as annotation strips, or a data.frame
#'   with an `observation_id` column plus the columns to display.
#' @param zscore Rescale each program's scores to mean 0 / SD 1 before
#'   plotting (default TRUE; keeps the color scale comparable across programs).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' plot_program_scores(fit)
#' # with an annotation strip built from a planted program's donor scores
#' sev <- sim$truth$B[[1]]$z > median(sim$truth$B[[1]]$z)
#' md <- data.frame(observation_id = paste0("obs", seq_along(sev)),
#'                  severity = factor(sev, c(FALSE, TRUE), c("low", "high")))
#' fit2 <- as_cell_program_data(sim$pseudobulk, sample_metadata = md)
#' plot_program_scores(canonicalize_programs(fit_celltype_programs(fit2, max_factors = 3, seed = 1)),
#'                     annotation = "severity")
#' @export
plot_program_scores <- function(fit, cell_types = NULL, annotation = NULL, zscore = TRUE) {
  cts <- (cell_types %||% names(fit$scores))
  cts <- cts[vapply(cts, function(ct) ncol(fit$scores[[ct]]) > 0L, TRUE)]
  if (!length(cts)) return(.empty_plot("No programs retained: nothing to plot."))
  S <- do.call(cbind, lapply(cts, function(ct) fit$scores[[ct]]))
  if (zscore) {
    sds <- apply(S, 2, stats::sd); sds[sds == 0] <- 1
    S <- scale(S, center = TRUE, scale = sds)
  }
  lims <- max(abs(S), na.rm = TRUE); if (!is.finite(lims) || lims == 0) lims <- 1
  df <- .melt_matrix(S)
  names(df) <- c("donor", "program", "score")
  df$donor <- factor(df$donor, levels = rownames(S))
  ylv <- rev(colnames(S))
  df$program <- factor(df$program, levels = ylv)
  df$score <- .clip_sym(df$score, lims)
  df$panel <- "program scores"

  ann <- NULL
  if (!is.null(annotation)) {
    if (is.character(annotation)) {
      md <- fit$sample_metadata[match(rownames(S), fit$sample_metadata$observation_id), , drop = FALSE]
      ann <- md[, annotation, drop = FALSE]
    } else {
      i <- match(rownames(S), annotation$observation_id)
      ann <- annotation[i, setdiff(names(annotation), "observation_id"), drop = FALSE]
    }
    rownames(ann) <- rownames(S)
    ann_cols <- setdiff(names(ann), "donor")
    long <- do.call(rbind, lapply(ann_cols, function(cl) {
      data.frame(donor = rownames(S), program = cl, score = NA_real_,
                 panel = cl, value = as.character(ann[, cl]), stringsAsFactors = FALSE)
    }))
    long$donor <- factor(long$donor, levels = rownames(S))
    ## annotation values are categorical: draw as fixed-fill tiles with text,
    ## on a fill scale separate from the numeric score gradient
    ylv <- c(ylv, ann_cols)
    df$program <- factor(df$program, levels = ylv)
    long$program <- factor(long$program, levels = ylv)
    df$value <- NA_character_
    df <- rbind(df[, c("donor", "program", "score", "panel", "value")], long)
    df$panel <- factor(df$panel, levels = c("program scores", ann_cols))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(donor, program, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = .diverging_pal()[1], mid = .diverging_pal()[2],
                                  high = .diverging_pal()[3], midpoint = 0, na.value = "grey95",
                                  limits = c(-lims, lims),
                                  name = if (zscore) "z-score" else "score") +
    ggplot2::labs(x = NULL, y = NULL) +
    .theme_programs() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  if (is.null(annotation)) {
    p
  } else {
    ann_df <- df[df$panel != "program scores", , drop = FALSE]
    p +
      ggplot2::geom_tile(data = ann_df, ggplot2::aes(donor, program),
                         fill = "grey92", color = "white", linewidth = 0.2) +
      ggplot2::geom_text(data = ann_df[ann_df$value != "NA" & !is.na(ann_df$value), , drop = FALSE],
                         ggplot2::aes(donor, program, label = value), size = 2.4, angle = 90) +
      ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y")
  }
}

## ---------------------------------------------------------------------------
## 2. top genes per program
## ---------------------------------------------------------------------------

#' Plot top-loaded genes per program
#'
#' Bar plot of the largest-magnitude loadings for one or more programs,
#' colored by loading sign. Loadings are posterior means: small nonzero
#' values are shrinkage noise, so extract gene sets with an explicit
#' threshold or lfsr rather than exact zeros.
#'
#' @param fit A canonicalized `cell_program_fit` object.
#' @param programs Program id(s) (e.g. `"B_1"` or `c("B_1", "CD4_1")`).
#' @param n Number of top genes per program (by |loading|).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' plot_top_genes(fit, c("B_1", "CD4_1"), n = 8)
#' @export
plot_top_genes <- function(fit, programs, n = 10) {
  rows <- lapply(programs, function(p) {
    ct <- sub("_[^_]+$", "", p)
    W <- fit$loadings[[ct]]
    if (is.null(W) || !p %in% colnames(W)) return(NULL)
    w <- W[, p]
    k <- order(-abs(w))[seq_len(min(n, length(w)))]
    data.frame(program = p, gene = names(w[k]), loading = unname(w[k]),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || !nrow(df)) return(.empty_plot("No matching programs with loadings."))
  df$gene <- factor(df$gene, levels = unique(df$gene[order(df$program, df$loading)]))
  df$program <- factor(df$program, levels = programs)
  ggplot2::ggplot(df, ggplot2::aes(loading, gene, fill = loading > 0)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#2166AC", `TRUE` = "#B2182B"),
                               labels = c(`FALSE` = "negative", `TRUE` = "positive"), name = NULL) +
    ggplot2::facet_wrap(~program, scales = "free_y") +
    ggplot2::labs(x = "loading", y = NULL) +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 3. variance explained
## ---------------------------------------------------------------------------

#' Plot variance explained per program
#'
#' Tile plot of the fraction of variance each program explains in its own
#' cell type's pseudobulk matrix (from the underlying flashier fit).
#'
#' @param fit A `cell_program_fit` object.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' plot_variance_explained(fit)
#' @export
plot_variance_explained <- function(fit) {
  rows <- lapply(names(fit$fits), function(ct) {
    fl <- fit$fits[[ct]]$flash
    pve <- if (!is.null(fl) && !is.null(fl$pve)) fl$pve else return(NULL)
    kn <- ncol(fit$scores[[ct]])
    if (!length(pve) || kn == 0L) return(NULL)
    m <- min(kn, length(pve))
    data.frame(cell_type = ct, program = colnames(fit$scores[[ct]])[seq_len(m)],
               pve = pve[seq_len(m)], stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || !nrow(df)) return(.empty_plot("No fitted programs: variance explained unavailable."))
  df$program <- factor(df$program, levels = unique(df$program))
  ggplot2::ggplot(df, ggplot2::aes(cell_type, program, fill = pve)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * pve)), size = 3) +
    ggplot2::scale_fill_viridis_c(name = "PVE", option = "mako", direction = -1) +
    ggplot2::labs(x = NULL, y = NULL) +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 4. sharing heatmap
## ---------------------------------------------------------------------------

#' Plot program sharing results
#'
#' Heatmap of pairwise cross-cell-type program correlations (the correlation
#' pre-screen in `program_sharing`), with each program's final mashr-based
#' classification appended to the axis labels.
#'
#' @param share Output of `program_sharing`.
#' @param level `"program"` (per-program matrix; readable up to a few dozen
#'   programs) or `"cell_type"` (mean pre-screen correlation per cell-type
#'   pair; the cohort-scale view).
#' @param q_alpha Significance level for marking a match (default 0.05).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' if (requireNamespace("mashr", quietly = TRUE)) { # program_sharing needs mashr (Suggests)
#'   share <- program_sharing(fit, n_perm = 200, seed = 1)
#'   plot_sharing(share)
#' }
#' @export
plot_sharing <- function(share, level = c("program", "cell_type"), q_alpha = 0.05) {
  level <- match.arg(level)
  m <- share$matches
  summ <- share$summary
  if (is.null(m) || !nrow(m)) {
    return(.empty_plot("No cross-cell-type program matches passed the pre-screen."))
  }
  if (level == "cell_type") {
    ## cohort-scale view: mean pre-screen correlation per cell-type pair
    agg <- stats::aggregate(cor ~ ct_a + ct_b, data = m,
                            FUN = function(x) max(abs(x)))
    both <- rbind(agg, stats::setNames(agg, c("ct_b", "ct_a", "cor")))
    cts <- sort(unique(c(m$ct_a, m$ct_b)))
    mat <- matrix(NA_real_, length(cts), length(cts), dimnames = list(cts, cts))
    mat[cbind(both$ct_a, both$ct_b)] <- both$cor
    df <- .melt_matrix(mat); names(df) <- c("a", "b", "cor")
    df$a <- factor(df$a, levels = cts)
    df$b <- factor(df$b, levels = rev(cts))
    return(ggplot2::ggplot(df, ggplot2::aes(a, b, fill = cor)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(cor), "", sprintf("%.2f", cor))),
                         size = 3) +
      ggplot2::scale_fill_gradient2(low = .diverging_pal()[1], mid = .diverging_pal()[2],
                                    high = .diverging_pal()[3], midpoint = 0,
                                    na.value = "grey95", name = "max |cor|") +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = NULL,
                    subtitle = "strongest cross-cell-type pre-screen correlation per cell-type pair") +
      .theme_programs())
  }
  progs <- unique(c(m$prog_a, m$prog_b))
  mat <- matrix(NA_real_, length(progs), length(progs), dimnames = list(progs, progs))
  qmat <- mat
  for (i in seq_len(nrow(m))) {
    mat[m$prog_a[i], m$prog_b[i]] <- m$cor[i]
    mat[m$prog_b[i], m$prog_a[i]] <- m$cor[i]
    qmat[m$prog_a[i], m$prog_b[i]] <- m$q[i]
    qmat[m$prog_b[i], m$prog_a[i]] <- m$q[i]
  }
  diag(mat) <- 1
  df <- .melt_matrix(mat); names(df) <- c("a", "b", "cor")
  df$q <- qmat[cbind(match(df$a, progs), match(df$b, progs))]
  ## asterisk labels swamp the plot beyond a few hundred matches
  df$lab <- if (nrow(m) <= 500L) ifelse(is.na(df$q), "", ifelse(df$q < q_alpha, "*", "")) else ""
  df$cor[is.na(df$cor) & df$a != df$b] <- NA_real_

  cls <- summ$classification[match(progs, summ$program_id)]
  cls[is.na(cls)] <- "not tested"
  ylab <- sprintf("%s (%s)", progs, cls)
  names(ylab) <- progs
  df$a <- factor(df$a, levels = progs, labels = ylab)
  df$b <- factor(df$b, levels = rev(progs), labels = rev(ylab))

  ggplot2::ggplot(df, ggplot2::aes(a, b, fill = cor)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 5, vjust = 0.8) +
    ggplot2::scale_fill_gradient2(low = .diverging_pal()[1], mid = .diverging_pal()[2],
                                  high = .diverging_pal()[3], midpoint = 0,
                                  na.value = "grey95", name = "cor") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL,
                  subtitle = "* q < 0.05 in the correlation pre-screen; axis labels show mashr classification") +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 5. sharing spectrum
## ---------------------------------------------------------------------------

#' Plot the sharing spectrum
#'
#' Per-pair summary of `sharing_spectrum`: the smallest principal angle
#' between each cell-type pair's program-activity subspaces, point size =
#' number of shared subspace dimensions (cosine > 0.5), color = permutation
#' test outcome. Hollow points are underpowered pairs.
#'
#' @param spectrum Output of `sharing_spectrum`.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' plot_sharing_spectrum(sharing_spectrum(fit, n_perm = 200))
#' @export
plot_sharing_spectrum <- function(spectrum) {
  s <- spectrum$summary
  if (is.null(s) || !nrow(s)) return(.empty_plot("Empty sharing spectrum."))
  s$label <- ifelse(s$underpowered, "underpowered",
                    ifelse(isTRUE(s$shared), sprintf("q = %.3f", s$q_pair), "n.s."))
  ggplot2::ggplot(s, ggplot2::aes(min_angle_deg, pair)) +
    ggplot2::geom_point(ggplot2::aes(size = n_shared_05, color = label, shape = underpowered),
                        stroke = 1.2) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
    ggplot2::scale_color_brewer(palette = "Set1", name = NULL) +
    ggplot2::scale_size_continuous(name = "shared dims\n(cos > 0.5)", range = c(3, 8)) +
    ggplot2::labs(x = "smallest principal angle (deg)", y = NULL) +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 6. artifact flags
## ---------------------------------------------------------------------------

#' Plot artifact-gene flags
#'
#' Bar plot of the fraction of each program's top genes that are V(D)J or
#' donor-genotype-marker genes, with the flag threshold marked.
#'
#' @param flags Output of `flag_artifact_genes`.
#' @param threshold Flag threshold (drawn as a dashed line; default 0.2).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' plot_artifact_flags(flag_artifact_genes(fit))
#' @export
plot_artifact_flags <- function(flags, threshold = 0.2) {
  if (is.null(flags) || !nrow(flags)) return(.empty_plot("No programs to flag."))
  df <- data.frame(program = paste(flags$cell_type, flags$program_id, sep = "_"),
                   vdj = flags$vdj_frac, genotype = flags$genotype_marker_frac,
                   vdj_flag = flags$vdj_flag, genotype_flag = flags$genotype_marker_flag,
                   stringsAsFactors = FALSE)
  long <- do.call(rbind, lapply(seq_len(nrow(df)), function(i)
    data.frame(program = df$program[i],
               set = c("V(D)J", "genotype marker"),
               fraction = c(df$vdj[i], df$genotype[i]),
               flagged = c(df$vdj_flag[i], df$genotype_flag[i]),
               stringsAsFactors = FALSE)))
  ggplot2::ggplot(long, ggplot2::aes(fraction, program, fill = set)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_brewer(palette = "Set1", name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "fraction of top genes", y = NULL,
                  subtitle = "dashed line: flag threshold") +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 7. CMP network
## ---------------------------------------------------------------------------

#' Plot the coordinated-program network
#'
#' Network of significant cross-cell-type program edges. Node size = node
#' statistic A_v, node color = CMP membership, edge width = |T|, edge color =
#' edge sign. Programs flagged as ambient get hollow nodes.
#'
#' @param cmp Output of `find_coordinated_programs`.
#' @param max_edges Draw only the strongest edges (default: all).
#' @param seed Layout seed (the spring layout is stochastic).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
#' plot_cmp_network(cmp)
#' @export
plot_cmp_network <- function(cmp, max_edges = NULL, seed = 1) {
  e <- cmp$edges; n <- cmp$nodes
  if (is.null(e) || !nrow(e)) return(.empty_plot("No significant coordination edges."))
  if (!is.null(max_edges) && nrow(e) > max_edges) {
    e <- e[order(-abs(e$stat)), ][seq_len(max_edges), , drop = FALSE]
  }
  memb <- stats::setNames(as.character(cmp$cmps$cmp), cmp$cmps$program_id)
  amb <- stats::setNames(
    n$ambient_score >= (cmp$settings$ambient_threshold %||% 0.3),
    paste(n$cell_type, n$program_id, sep = "_"))

  nodes <- unique(c(e$prog_a, e$prog_b))
  nd <- data.frame(label = nodes, cmp_lab = ifelse(is.na(memb[nodes]), "unassigned",
                                                   paste("CMP", memb[nodes])),
                   ambient = unname(amb[nodes]) %in% TRUE,
                   node_stat = n$node_stat[match(nodes, n$program_id)],
                   stringsAsFactors = FALSE)
  nd$node_stat[is.na(nd$node_stat)] <- 0

  ed <- data.frame(from = e$prog_a, to = e$prog_b,
                   stat = e$stat, sign = e$sign, stringsAsFactors = FALSE)

  pos <- tryCatch({
    if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph not installed")
    set.seed(seed)
    g <- igraph::graph_from_data_frame(ed[, c("from", "to")], directed = FALSE,
                                       vertices = data.frame(name = nd$label))
    lay <- igraph::layout_with_fr(g)
    data.frame(label = nd$label, x = lay[, 1], y = lay[, 2], stringsAsFactors = FALSE)
  }, error = function(err) {
    th <- seq(0, 2 * pi, length.out = length(nodes) + 1)[-1]
    data.frame(label = nodes, x = cos(th), y = sin(th), stringsAsFactors = FALSE)
  })
  nd <- merge(nd, pos, by = "label")
  ed <- merge(ed, nd[, c("label", "x", "y")], by.x = "from", by.y = "label")
  names(ed)[names(ed) %in% c("x", "y")] <- c("x0", "y0")
  ed <- merge(ed, nd[, c("label", "x", "y")], by.x = "to", by.y = "label")
  names(ed)[names(ed) %in% c("x", "y")] <- c("x1", "y1")

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed,
                          ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1,
                                       linewidth = abs(stat), color = factor(sign)),
                          lineend = "round") +
    ggplot2::geom_point(data = nd,
                        ggplot2::aes(x, y, size = node_stat, fill = cmp_lab, shape = ambient),
                        color = "black", stroke = 0.8) +
    ggplot2::geom_text(data = nd, ggplot2::aes(x, y, label = label), size = 2.8,
                       nudge_y = 0.12, check_overlap = TRUE) +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 2), name = "|T|") +
    ggplot2::scale_color_manual(values = c(`1` = "#B2182B", `-1` = "#2166AC"),
                                name = "edge sign", labels = c(`1` = "+", `-1` = "-")) +
    ggplot2::scale_size_continuous(name = "node stat", range = c(4, 10)) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
    ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    .theme_programs() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
}

## ---------------------------------------------------------------------------
## 8. CMP node statistics
## ---------------------------------------------------------------------------

#' Plot CMP node statistics
#'
#' Bar plot of each program's node statistic A_v (pooled coordination
#' evidence across partner cell types), colored by CMP membership, with
#' permutation p-values marked.
#'
#' @param cmp Output of `find_coordinated_programs`.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
#' plot_cmp_nodes(cmp)
#' @export
plot_cmp_nodes <- function(cmp) {
  n <- cmp$nodes
  if (is.null(n) || !nrow(n)) return(.empty_plot("No programs entered the coordination test."))
  n$label <- paste(n$cell_type, n$program_id, sep = "_")
  memb <- stats::setNames(as.character(cmp$cmps$cmp), cmp$cmps$program_id)
  n$cmp_lab <- ifelse(is.na(memb[n$program_id]), "unassigned", paste("CMP", memb[n$program_id]))
  n$sig <- ifelse(n$node_p < (cmp$settings$alpha %||% 0.05), "*", "")
  ggplot2::ggplot(n, ggplot2::aes(node_stat, stats::reorder(label, node_stat))) +
    ggplot2::geom_col(ggplot2::aes(fill = cmp_lab), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sig), vjust = 0.4, hjust = -0.3, size = 5) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = "node statistic (pooled |T| excess)", y = NULL,
                  subtitle = "* permutation p < 0.05") +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 8b. CMP membership across settings (alluvial)
## ---------------------------------------------------------------------------

#' Plot CMP membership across settings
#'
#' Alluvial view of how each program's coordination membership changes across
#' a named list of `cmp_result` objects — typically the same fit run with
#' different `find_coordinated_programs()` settings, e.g. `ambient = "drop"`
#' versus `ambient = "flag"`. Strata are membership classes ("CMP k", 
#' "significant, no CMP", "not significant", "excluded (ambient)"); ribbons
#' are programs, colored by their destination class, so a setting change that
#' calls or drops programs shows up as a colored ribbon leaving the grey
#' baseline.
#'
#' @param cmp_list A named list of `cmp_result` objects; the names become the
#'   x-axis labels, in order.
#' @param alpha Membership significance level used to classify programs
#'   (default 0.05; each result's own `settings$alpha` takes precedence).
#' @param label_programs When TRUE and at most 12 programs are present,
#'   annotate ribbons whose membership changed with the program IDs carried.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' cmps <- list(
#'   drop = find_coordinated_programs(fit, n_perm = 200, seed = 1),
#'   flag = find_coordinated_programs(fit, n_perm = 200, seed = 1, ambient = "flag")
#' )
#' plot_cmp_membership(cmps)
#' @export
plot_cmp_membership <- function(cmp_list, alpha = 0.05, label_programs = FALSE) {
  if (!is.list(cmp_list) || is.null(names(cmp_list)) ||
      any(!nzchar(names(cmp_list))) ||
      !all(vapply(cmp_list, inherits, TRUE, "cmp_result"))) {
    .stopf("`cmp_list` must be a named list of cmp_result objects.")
  }
  n_s <- length(cmp_list)

  ## -- classify every program of one result
  classify <- function(cmp) {
    a <- if (!is.null(cmp$settings$alpha)) cmp$settings$alpha else alpha
    out <- character(0)
    if (!is.null(cmp$nodes) && nrow(cmp$nodes)) {
      k <- paste(cmp$nodes$cell_type, cmp$nodes$program_id, sep = "|")
      memb <- stats::setNames(as.character(cmp$cmps$cmp), cmp$cmps$program_id)
      m <- memb[cmp$nodes$program_id]
      lv <- ifelse(is.na(m), "significant, no CMP", paste("CMP", m))
      lv[cmp$nodes$node_p >= a] <- "not significant"
      out[k] <- lv
    }
    if (!is.null(cmp$ambient) && nrow(cmp$ambient)) {
      k <- paste(cmp$ambient$cell_type, cmp$ambient$program_id, sep = "|")
      thr <- if (!is.null(cmp$settings$ambient_threshold)) cmp$settings$ambient_threshold else 0.3
      exc <- !k %in% names(out) & cmp$ambient$ambient_score > thr
      out[k[exc]] <- "excluded (ambient)"
      out[setdiff(k, names(out))] <- "not tested"
    }
    out
  }
  cls <- lapply(cmp_list, classify)

  ## -- canonicalize CMP labels across settings by member-set overlap: the
  ##    CMP index from find_coordinated_programs is arbitrary per run, so a
  ##    CMP that keeps its members should keep its label
  prev_members <- list()
  for (i in seq_len(n_s)) {
    lv <- cls[[i]]
    is_cmp <- grepl("^CMP ", lv)
    raw <- unique(lv[is_cmp])
    raw <- raw[order(as.integer(sub("CMP ", "", raw)))]
    mapping <- character(0); used <- character(0)
    for (k in raw) {
      mem <- names(lv)[lv == k]
      best_j <- 0; best <- NA_character_
      for (lab in names(prev_members)) {
        if (lab %in% used) next
        j <- length(intersect(mem, prev_members[[lab]])) /
          max(length(union(mem, prev_members[[lab]])), 1L)
        if (j > best_j) { best_j <- j; best <- lab }
      }
      if (best_j > 0) {
        mapping[k] <- best
      } else {
        used_nums <- as.integer(sub("CMP ", "", c(names(prev_members), used)))
        mapping[k] <- paste("CMP", min(setdiff(seq_len(max(c(used_nums, 0L)) + 1L),
                                               used_nums)))
      }
      used <- c(used, mapping[k])
    }
    lv[is_cmp] <- mapping[lv[is_cmp]]
    cls[[i]] <- lv
    prev_members <- setNames(lapply(mapping, function(k) names(lv)[lv == k]),
                             mapping)
  }

  keys <- sort(unique(unlist(lapply(cls, names))))
  mat <- matrix("not tested", length(keys), n_s,
                dimnames = list(keys, names(cmp_list)))
  for (i in seq_len(n_s)) mat[names(cls[[i]]), i] <- cls[[i]]

  cmp_n <- max(as.integer(sub("CMP ", "", grep("^CMP ", mat, value = TRUE))), 0L)
  lvls <- c(if (cmp_n > 0) paste("CMP", seq_len(cmp_n)), "significant, no CMP",
            "not significant", "excluded (ambient)", "not tested")

  ## -- strata rectangles per setting (stacked bottom-up in level order)
  bw <- 0.14
  rects <- do.call(rbind, lapply(seq_len(n_s), function(i) {
    tab <- table(factor(mat[, i], levels = lvls))
    tab <- tab[tab > 0]
    y0 <- 0
    do.call(rbind, lapply(names(tab), function(lv) {
      y1 <- y0 + tab[[lv]]
      y0 <<- y1
      data.frame(x = i, ymin = y1 - tab[[lv]], ymax = y1, level = lv,
                 n = as.integer(tab[[lv]]), stringsAsFactors = FALSE)
    }))
  }))
  y_start <- function(i, lv) {
    r <- rects[rects$x == i & rects$level == lv, ]
    if (!nrow(r)) NA_real_ else r$ymin[1]
  }

  ## -- flow rectangles between consecutive settings, allocated in level order
  flows <- do.call(rbind, lapply(seq_len(n_s - 1), function(i) {
    L <- mat[, i]; R <- mat[, i + 1]
    off_L <- off_R <- numeric(0)
    do.call(rbind, lapply(lvls, function(lf) {
      do.call(rbind, lapply(lvls, function(lt) {
        ks <- names(L)[L == lf & R == lt]
        n <- length(ks)
        if (!n) return(NULL)
        a0 <- y_start(i, lf); b0 <- y_start(i + 1, lt)
        if (is.na(a0) || is.na(b0)) return(NULL)
        a1 <- if (is.na(off_L[lf])) a0 else off_L[lf]
        b1 <- if (is.na(off_R[lt])) b0 else off_R[lt]
        off_L[lf] <<- a1 + n; off_R[lt] <<- b1 + n
        data.frame(xmin = i + bw, xmax = i + 1 - bw, ymin = a1, ymax = a1 + n,
                   ymin2 = b1, ymax2 = b1, from = lf, to = lt,
                   programs = paste(sub("^[^|]*\\|", "", ks), collapse = ", "),
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

  ## -- colors: Dark2 for CMPs, greys for the rest
  has_flows <- !is.null(flows) && nrow(flows) > 0L
  lvls_present <- lvls[lvls %in% unique(c(rects$level, if (has_flows) flows$to))]
  n_cmp <- sum(grepl("^CMP ", lvls_present))
  pal <- c(setNames(grDevices::hcl.colors(max(n_cmp, 1L), "Dark2"),
                    grep("^CMP ", lvls_present, value = TRUE)),
           setNames(c("grey72", "grey40", "grey88"),
                    c("significant, no CMP", "excluded (ambient)", "not tested")))
  pal <- pal[lvls_present]

  p <- ggplot2::ggplot()
  if (has_flows) {
    p <- p +
      ggplot2::geom_rect(data = flows, ggplot2::aes(xmin = xmin, xmax = xmax,
                                                    ymin = ymin, ymax = ymax,
                                                    fill = to), color = NA,
                         alpha = 0.55) +
      ggplot2::geom_rect(data = flows, ggplot2::aes(xmin = xmin, xmax = xmax,
                                                    ymin = ymin2, ymax = ymax2,
                                                    fill = to), color = NA,
                         alpha = 0.55)
  }
  p <- p +
    ggplot2::geom_rect(data = rects, ggplot2::aes(xmin = x - bw, xmax = x + bw,
                                                  ymin = ymin, ymax = ymax,
                                                  fill = level),
                       color = "white", linewidth = 0.3) +
    ggplot2::geom_text(data = rects, ggplot2::aes(x = x, y = (ymin + ymax) / 2,
                                                  label = paste0(level, " (", n, ")"),
                                                  hjust = ifelse(x == 1, 1.12, -0.12)),
                       size = 3.1) +
    ggplot2::scale_fill_manual(values = pal, breaks = lvls_present, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_s), labels = names(cmp_list),
                                expand = ggplot2::expansion(mult = 0.35)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(x = NULL, y = "programs") +
    .theme_programs() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                   legend.position = "none",
                   panel.grid = ggplot2::element_blank())
  if (label_programs && length(keys) <= 12L && has_flows) {
    ch <- flows[flows$from != flows$to, ]
    if (nrow(ch)) {
      p <- p + ggplot2::geom_text(data = ch, ggplot2::aes(x = (xmin + xmax) / 2,
                                                          y = (ymin2 + ymax2) / 2,
                                                          label = programs),
                                  size = 2.8, hjust = 0.5, vjust = 0.5)
    }
  }
  p
}

## ---------------------------------------------------------------------------
## 9. context association
## ---------------------------------------------------------------------------

#' Plot CMP context associations
#'
#' Dot-and-whisker plot of the context-variable effect on each CMP's donor
#' scores (linear model from `cmp_context_association`), with 95% CIs from
#' the model standard errors.
#'
#' @param assoc Output of `cmp_context_association` (data.frame with `cmp`,
#'   `term`, `estimate`, `t`, `p`, `n`).
#' @param p_alpha Significance level for solid vs open points (default 0.05).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
#' # context variable: high/low donor score of the first planted program
#' z <- sim$truth$B[[1]]$z > median(sim$truth$B[[1]]$z)
#' ctx <- setNames(factor(z, c(FALSE, TRUE), c("low", "high")), paste0("obs", seq_along(z)))
#' plot_context_association(cmp_context_association(fit, cmp, ctx))
#' @export
plot_context_association <- function(assoc, p_alpha = 0.05) {
  if (is.null(assoc) || !nrow(assoc)) return(.empty_plot("No context associations to plot."))
  df <- assoc
  df$se <- abs(df$estimate / df$t)
  df$cmp_lab <- paste0("CMP", df$cmp)
  ## the context term is named "x" by the model formula; label it plainly
  df$term <- ifelse(df$term == "x", "context variable", df$term)
  ggplot2::ggplot(df, ggplot2::aes(estimate, cmp_lab, color = p < p_alpha)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                           orientation = "y", width = 0.15, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = p < p_alpha), size = 3) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey40"),
                                name = sprintf("p < %.2f", p_alpha)) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    ggplot2::facet_wrap(~term, scales = "free_x") +
    ggplot2::labs(x = "effect on CMP score (per context unit)", y = NULL) +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 10. alignment / principal angles
## ---------------------------------------------------------------------------

#' Plot program alignment (principal angles)
#'
#' Visualize `principal_angles` / `align_programs` results. For a single
#' `principal_angles` object: how each program of each cell type loads on the
#' shared principal directions (program-weight tiles), annotated with the
#' permutation-based subspace test. For an `align_programs` result (named
#' list of such objects): faceted per cell-type pair.
#'
#' @param x A `principal_angles` object, or the named list returned by
#'   `align_programs`.
#' @param q_alpha Significance level for the per-pair subspace test
#'   (default 0.05).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' plot_alignment(principal_angles(fit, "B", "CD4", n_perm = 200))
#' @export
plot_alignment <- function(x, q_alpha = 0.05) {
  if (inherits(x, "principal_angles")) return(.plot_one_alignment(x, q_alpha))
  ref <- attr(x, "reference")
  rows <- lapply(names(x), function(ct) {
    r <- x[[ct]]
    if (!inherits(r, "principal_angles")) return(NULL)
    df <- .alignment_df(r)
    df$pair <- paste(ref, "vs", ct)
    df
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || !nrow(df)) return(.empty_plot("No alignment results to plot."))
  sig_lab <- if (any(df$q_pair < q_alpha, na.rm = TRUE)) {
    sprintf("%d/%d pairs significant (q < %.2f)", sum(df$q_pair < q_alpha, na.rm = TRUE),
            length(unique(df$pair)), q_alpha)
  } else "no pair significant"
  ggplot2::ggplot(df, ggplot2::aes(direction, program, fill = weight)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = .diverging_pal()[1], mid = .diverging_pal()[2],
                                  high = .diverging_pal()[3], midpoint = 0,
                                  name = "principal-\nvector weight") +
    ggplot2::facet_grid(pair ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "shared principal direction", y = NULL, subtitle = sig_lab) +
    .theme_programs()
}

.alignment_df <- function(r) {
  pa <- r$program_weights$a; pb <- r$program_weights$b
  mk <- function(M) {
    if (is.null(M) || !nrow(M)) return(NULL)
    d <- .melt_matrix(M); names(d) <- c("program", "direction", "weight")
    ## both sides must share one direction index (their column names differ)
    d$direction <- match(d$direction, colnames(M))
    d
  }
  df <- rbind(mk(pa), mk(pb))
  if (is.null(df)) return(df)
  df$direction <- factor(df$direction, levels = sort(unique(df$direction)),
                         labels = paste("dir", sort(unique(df$direction))))
  df$q_pair <- r$q_pair
  df
}

.plot_one_alignment <- function(r, q_alpha) {
  df <- .alignment_df(r)
  if (is.null(df) || !nrow(df)) return(.empty_plot("No shared principal directions."))
  ## a bare principal_angles object carries the raw permutation p (q is added
  ## by sharing_spectrum's BH step); fall back through q -> p
  qp <- if (!is.null(r$q_pair) && length(r$q_pair) == 1L && !is.na(r$q_pair)) r$q_pair else r$p_pair
  sig <- if (isTRUE(r$underpowered) || is.null(qp) || is.na(qp)) "underpowered" else
    if (qp < q_alpha) sprintf("p = %.3f", qp) else "n.s."
  df$program <- factor(df$program, levels = rev(unique(df$program)))
  ggplot2::ggplot(df, ggplot2::aes(direction, program, fill = weight)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = .diverging_pal()[1], mid = .diverging_pal()[2],
                                  high = .diverging_pal()[3], midpoint = 0,
                                  name = "principal-\nvector weight") +
    ggplot2::labs(x = "shared principal direction", y = NULL,
                  subtitle = paste("subspace test:", sig)) +
    .theme_programs()
}

## ---------------------------------------------------------------------------
## 11. stability
## ---------------------------------------------------------------------------

#' Plot program stability
#'
#' Dot plot of bootstrap stability per program: recovery frequency (fraction
#' of bootstrap refits in which the program reappears) against the median
#' score correlation between the original and matched bootstrap programs.
#'
#' @param fit A `cell_program_fit` object processed by
#'   `assess_program_stability` (uses its `$stability` element).
#' @return A ggplot object.
#' @examples
#' sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
#' x <- as_cell_program_data(sim$pseudobulk)
#' fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
#' fit <- assess_program_stability(fit, n_boot = 20, seed = 1)
#' plot_stability(fit)
#' @export
plot_stability <- function(fit) {
  st <- fit$stability$per_program
  if (is.null(st) || !nrow(st)) {
    return(.empty_plot("No stability results: run assess_program_stability() first."))
  }
  ggplot2::ggplot(st, ggplot2::aes(recovery_freq, median_score_cor)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(size = top_gene_jaccard, color = sign_consistency),
                        alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = .diverging_pal()[1], mid = "grey90",
                                   high = .diverging_pal()[3], midpoint = 0,
                                   limits = c(-1, 1), name = "sign\nconsistency") +
    ggplot2::scale_size_continuous(name = "top-gene\nJaccard", range = c(2, 8)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "bootstrap recovery frequency", y = "median score correlation",
                  title = "Program stability across bootstrap refits") +
    .theme_programs()
}
