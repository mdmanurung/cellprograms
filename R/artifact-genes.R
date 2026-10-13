# cellprograms: artifact-gene sets, CMP flagging, and fit-time exclusion
#
# Real-data CMP detection at max_factors = 30 roughly doubles CMP counts, but
# many of the additional CMPs are carried by V(D)J genes or by genes with
# common germline variation (deletions, copy-number polymorphisms, HLA
# haplotypes, sex chromosomes). Such genes co-vary across donors because of
# genotype or clonal composition, not coordinated biology. This file provides
# (a) the gene sets (`artifact_gene_sets`), (b) flagging of called CMPs or
# programs that carry them (`flag_artifact_genes()`), and (c) the universe
# expansion used for fit-time exclusion (`artifact_genes()`; see the
# `exclude_genes` argument of `fit_celltype_programs()`).

#' Artifact-gene sets
#'
#' Gene classes that can produce cross-donor co-variation without coordinated
#' biology, used for flagging called programs and for fit-time exclusion:
#'
#' * `vdj`: immunoglobulin and T-cell-receptor V(D)J and constant genes.
#'   Stored as a single regular expression,
#'   `^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)`, not an enumerated gene
#'   list, because the full V-gene complement is genome-annotation-dependent;
#'   expand it against a gene universe with `.expand_vdj()` (or use
#'   `artifact_genes()`, which expands every set).
#' * `genotype_marker`: genes with common germline deletion, copy-number or
#'   HLA-haplotype expression variation (GSTM1/GSTT1 deletions, CCL3L1/CCL4L1
#'   copy number, HLA-DRB3/4/5 and HLA-DQA/DQB haplotypes, ...).
#' * `sex`: XIST and the male-specific Y-chromosome genes RPS4Y1, DDX3Y, UTY,
#'   KDM5D, EIF1AY, NLGN4Y, TMSB4Y, USP9Y, ZFY.
#'
#' Provenance: the `genotype_marker` list is a screening heuristic, not
#' validated against genotype data; extend or replace it for your cohort. The
#' sets overlap `cmp_nuisance_genes()` (used by `cmp_triage()`); unlike that
#' regex vector, `artifact_gene_sets` stores explicit gene IDs where the set
#' is closed, so it can be subtracted from a feature panel directly.
#'
#' @format A named list of character vectors (`vdj` is a length-1 regular
#'   expression; `genotype_marker` and `sex` are gene IDs).
#' @export
artifact_gene_sets <- list(
  vdj = "^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)",
  genotype_marker = c("GSTM1", "GSTT1", "GSTT2", "UTS2", "ERAP2",
                      "CCL3L1", "CCL3L3", "CCL4L1", "CCL4L2",
                      "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
                      "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2"),
  sex = c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY",
          "NLGN4Y", "TMSB4Y", "USP9Y", "ZFY")
)

# Expand the V(D)J pattern against a gene universe. The set is stored as a
# regex (see artifact_gene_sets), so membership is annotation-dependent.
.expand_vdj <- function(genes) grep(artifact_gene_sets$vdj, genes, value = TRUE)

# Expand one artifact set against a gene universe: regex expansion for `vdj`,
# intersection for the enumerated sets.
.artifact_expand <- function(set, genes) {
  if (identical(set, "vdj")) return(.expand_vdj(genes))
  intersect(artifact_gene_sets[[set]], genes)
}

#' Expand artifact-gene sets against a gene universe
#'
#' Returns the union of the requested `artifact_gene_sets` expanded against
#' `genes` (the V(D)J regex is matched; enumerated sets are intersected). Use
#' it to build the `exclude_genes` argument of `fit_celltype_programs()`.
#'
#' @param genes Character vector of gene IDs (the universe, e.g. the column
#'   names of a pseudobulk matrix).
#' @param sets Subset of `names(artifact_gene_sets)` (default all).
#'
#' @return Character vector of gene IDs in `genes` belonging to any requested
#'   set, in universe order.
#'
#' @export
artifact_genes <- function(genes, sets = names(artifact_gene_sets)) {
  unknown <- setdiff(sets, names(artifact_gene_sets))
  if (length(unknown)) {
    cli::cli_abort("Unknown artifact gene set(s): {.val {unknown}}.")
  }
  hits <- unique(unlist(lapply(sets, .artifact_expand, genes = genes)))
  intersect(genes, hits)
}

#' Flag programs or CMPs carrying artifact genes
#'
#' For each CMP (or each program when `cmp` is `NULL`), takes the top-`top_n`
#' genes by absolute loading of each member program (pooled and de-duplicated
#' across members for a CMP) and reports, for each artifact set, the fraction
#' of those genes matching the set and a flag that is TRUE when at least one
#' gene matches. This mirrors the screening heuristic used in the project's
#' CMP report; a flag means "inspect before interpreting", not "discard".
#'
#' @param fit A canonicalized `cell_program_fit`.
#' @param cmp A `cmp_result` or its `cmps` table; `NULL` (default) flags every
#'   program of every cell type instead.
#' @param sets Subset of `names(artifact_gene_sets)` to screen against
#'   (default `c("vdj", "genotype_marker")`).
#' @param top_n Number of top-|loading| genes per program (default 12).
#'
#' @return Data frame. With `cmp`: one row per CMP (`cmp`, `n_members`).
#'   With `cmp = NULL`: one row per program (`cell_type`, `program_id`).
#'   Both have, per set `<s>`, columns `<s>_frac` (fraction of the top genes
#'   matching) and `<s>_flag` (TRUE when >= 1 gene matches).
#'
#' @export
flag_artifact_genes <- function(fit, cmp = NULL,
                                sets = c("vdj", "genotype_marker"),
                                top_n = 12L) {
  stopifnot(inherits(fit, "cell_program_fit"))
  unknown <- setdiff(sets, names(artifact_gene_sets))
  if (length(unknown)) {
    cli::cli_abort("Unknown artifact gene set(s): {.val {unknown}}.")
  }
  per_cmp <- !is.null(cmp)
  tab <- if (per_cmp) {
    .cmp_table(cmp)
  } else {
    do.call(rbind, lapply(fit$cell_types, function(ct) {
      progs <- colnames(fit$loadings[[ct]])
      if (length(progs) == 0L) return(NULL)
      data.frame(cmp = NA_integer_, cell_type = ct, program_id = progs,
                 stringsAsFactors = FALSE)
    }))
  }
  if (is.null(tab) || !nrow(tab)) {
    empty <- if (per_cmp) {
      data.frame(cmp = integer(0), n_members = integer(0), stringsAsFactors = FALSE)
    } else {
      data.frame(cell_type = character(0), program_id = character(0),
                 stringsAsFactors = FALSE)
    }
    for (s in sets) {
      empty[[paste0(s, "_frac")]] <- numeric(0)
      empty[[paste0(s, "_flag")]] <- logical(0)
    }
    return(empty)
  }

  top_genes <- function(ct, prog) {
    w <- fit$loadings[[ct]][, prog]
    names(sort(-abs(w)))[seq_len(min(top_n, length(w)))]
  }
  groups <- if (per_cmp) split(seq_len(nrow(tab)), tab$cmp) else
    stats::setNames(as.list(seq_len(nrow(tab))), seq_len(nrow(tab)))

  rows <- lapply(names(groups), function(k) {
    m <- tab[groups[[k]], , drop = FALSE]
    tops <- unique(unlist(Map(top_genes, m$cell_type, m$program_id)))
    row <- if (per_cmp) {
      data.frame(cmp = m$cmp[1], n_members = nrow(m), stringsAsFactors = FALSE)
    } else {
      data.frame(cell_type = m$cell_type[1], program_id = m$program_id[1],
                 stringsAsFactors = FALSE)
    }
    for (s in sets) {
      hits <- .artifact_expand(s, tops)
      row[[paste0(s, "_frac")]] <- length(hits) / max(length(tops), 1L)
      row[[paste0(s, "_flag")]] <- length(hits) >= 1L
    }
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
