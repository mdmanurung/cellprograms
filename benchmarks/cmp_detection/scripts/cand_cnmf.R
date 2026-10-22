## cand_cnmf.R -- MR-2 cNMF fit-level arm (PREREG_mr2-cnmf-v1.md).
## Replaces the Stage A flashier fit: per cell type, run python cNMF on the
## donor x gene counts (K grid {5,10,20,30}, registered stability rule), then
## emit the fit-object fields the caller consumes (scores = donor x K usage,
## loadings = gene x K spectra, non-negative; canonicalize_programs SKIPPED --
## no signed factors exist; recorded as a structural difference, not a bug).
## The caller (n_agg_noanc) is unchanged.

CNMF_PY <- Sys.getenv("CNMF_PY", "python3")
CNMF_SCRIPT <- Sys.getenv("CNMF_SCRIPT", local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  f <- if (length(a)) sub("^--file=", "", a)[1] else NA_character_
  file.path(if (!is.na(f)) dirname(normalizePath(f)) else getwd(), "cnmf_fit.py")
}))

cnmf_fit_object <- function(x, kgrid = c(5, 10, 20, 30), tmpdir = tempfile()) {
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  scores <- list(); loadings <- list(); ks <- list()
  for (ct in x$cell_types) {
    ## MR-2 input decision (registered deviation, see PREREG_mr2-cnmf-v1.md Addendum):
    ## the planted semi matrices are continuous (centered signal, min < 0), but cNMF
    ## requires non-negative integer counts. cNMF therefore consumes the planted
    ## donor x gene matrix rounded + clipped to >= 0 (its native count input), NOT
    ## the centered/residualized matrix K30 fits. This tests each backend on its
    ## natural input; the alternative (shift the processed matrix) fabricates an
    ## arbitrary offset and is methodologically worse.
    Y <- as.matrix(x$pseudobulk[[ct]])          # donors x genes, planted (continuous)
    Y <- pmax(round(Y), 0)                      # recover count space for cNMF
    storage.mode(Y) <- "integer"
    stopifnot(all(Y == round(Y)), min(Y) >= 0)
    pref <- file.path(tmpdir, paste0("ct_", gsub("[^A-Za-z0-9_]", "_", ct)))
    io <- .write_counts_mm(Y, pref)
    r <- system2(CNMF_PY, c(shQuote(CNMF_SCRIPT), shQuote(io), shQuote(pref),
                            "--kgrid", shQuote(paste(kgrid, collapse = ","))),
                 stdout = "|", stderr = "|")
    cat(r, "\n")
    usage <- read.delim(paste0(pref, ".usage.tsv"), row.names = 1, check.names = FALSE)
    spectra <- read.delim(paste0(pref, ".spectra.tsv"), row.names = 1, check.names = FALSE)
    ks[[ct]] <- as.integer(readLines(paste0(pref, ".chosen_k.txt")))
    ## align donor/gene ids to the data object
    usage <- usage[match(rownames(Y), rownames(usage)), , drop = FALSE]
    rownames(usage) <- rownames(Y)
    spectra <- spectra[match(colnames(Y), rownames(spectra)), , drop = FALSE]
    rownames(spectra) <- colnames(Y)
    colnames(usage) <- paste0("P", seq_len(ncol(usage)))
    colnames(spectra) <- paste0("P", seq_len(ncol(spectra)))
    scores[[ct]] <- as.matrix(usage)
    loadings[[ct]] <- as.matrix(spectra)
  }
  out <- list(cell_types = x$cell_types, observation_ids = x$observation_ids,
              scores = scores, loadings = loadings,
              canonicalization = lapply(x$cell_types, function(ct)
                list(scale_multiplier = 1, sign_multiplier = 1)) |>
                stats::setNames(x$cell_types),
              sample_metadata = x$sample_metadata,
              chosen_k = ks, backend = "cnmf")
  class(out) <- c("cell_program_fit", "cell_program_data", "list")
  out
}

.write_counts_mm <- function(Y, prefix) {
  io <- paste0(prefix, ".mtx")
  M <- Matrix::Matrix(Y, sparse = TRUE)
  Matrix::writeMM(M, io)          # writeMM appends ".mtx" -> "<prefix>.mtx.mtx"
  if (!file.exists(io) && file.exists(paste0(io, ".mtx")))
    file.rename(paste0(io, ".mtx"), io)
  io
}
