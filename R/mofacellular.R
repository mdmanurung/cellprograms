# cellprograms: export adapter for MOFAcellulaR
#
# Converts the canonical cell_program_data pseudobulk representation into the
# (counts, coldata) input pair expected by MOFAcellulaR::create_init_exp(),
# so a joint multicellular factor analysis (MOFA2, one view per cell type)
# can be run on exactly the same pseudobulk profiles that
# fit_celltype_programs() consumes locally.

#' Export pseudobulk profiles to MOFAcellulaR input format
#'
#' Builds the `counts` matrix (genes in rows, pseudobulk profiles in columns,
#' column names `<cell_type>_<observation_id>`) and the accompanying
#' `coldata` data.frame (`donor_id`, `cell_type`, `cell_counts`) expected by
#' `MOFAcellulaR::create_init_exp()`. Works on both `cell_program_data` and
#' `cell_program_fit` objects (the pseudobulk matrices are carried through
#' fitting), guaranteeing that local EBMF and joint MOFA2 analyses see
#' identical inputs.
#'
#' @param x A `cell_program_data` or `cell_program_fit` object.
#' @param cell_counts Optional data.frame with columns `observation_id`,
#'   `cell_type`, and `cell_counts` (number of cells behind each pseudobulk
#'   profile). If NULL (default), `cell_counts` is filled with NA; supply it
#'   if you intend to use MOFAcellulaR's cell-count-based QC
#'   (`filt_profiles(ncells = ...)`).
#'
#' @return A list with elements `counts` (genes x profiles matrix) and
#'   `coldata` (data.frame with rownames matching `colnames(counts)`),
#'   ready for `MOFAcellulaR::create_init_exp(counts = ..., coldata = ...)`.
#'
#' @details MOFAcellulaR's TMM normalization (`tmm_trns`) expects
#'   *raw-count* pseudobulk; pass raw-count pseudobulk in `x` if you plan to
#'   use it. Log-normalized pseudobulk is fine for the factor analysis itself
#'   but not for TMM.
#'
#' @export
as_mofacellular_input <- function(x, cell_counts = NULL) {
  stopifnot(inherits(x, "cell_program_data"))

  counts_list <- list()
  donor_list <- list()
  ct_list <- list()
  for (ct in x$cell_types) {
    m <- x$pseudobulk[[ct]]
    # genes x observations for this cell type
    counts_list[[ct]] <- t(m)
    donor_list[[ct]] <- rownames(m)
    ct_list[[ct]] <- rep(ct, nrow(m))
  }

  counts <- do.call(cbind, counts_list)
  colnames(counts) <- unlist(lapply(x$cell_types, function(ct) {
    paste0(ct, "_", donor_list[[ct]])
  }), use.names = FALSE)

  coldata <- data.frame(
    donor_id = unlist(donor_list, use.names = FALSE),
    cell_type = unlist(ct_list, use.names = FALSE),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  if (is.null(cell_counts)) {
    coldata$cell_counts <- NA_integer_
  } else {
    if (!is.data.frame(cell_counts) ||
        !all(c("observation_id", "cell_type", "cell_counts") %in% names(cell_counts))) {
      cli::cli_abort(
        "{.arg cell_counts} must be a data.frame with columns {.val observation_id}, {.val cell_type}, {.val cell_counts}."
      )
    }
    key <- paste0(coldata$cell_type, "_", coldata$donor_id)
    ckey <- paste0(cell_counts$cell_type, "_", cell_counts$observation_id)
    idx <- match(key, ckey)
    if (anyNA(idx)) {
      cli::cli_abort("{.arg cell_counts} is missing {sum(is.na(idx))} pseudobulk profile(s).")
    }
    coldata$cell_counts <- as.integer(cell_counts$cell_counts[idx])
  }

  list(counts = counts, coldata = coldata)
}
