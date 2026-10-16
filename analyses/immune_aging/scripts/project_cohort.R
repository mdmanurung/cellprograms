# Cross-cohort MR-1 projection + conservation.
#
# Projects a discovery-cohort fit (frozen preprocessing, canonicalized) onto a
# query cohort's donor x gene pseudobulks, per cell type (tryCatch-isolated),
# then reports:
#   - projection coverage per cell type (>= 0.9 gate)
#   - per-program concordance: cor(projected score, query refit score) across
#     shared donors, after optimal loading matching (query refit vs discovery)
#
# Usage: Rscript project_cohort.R <discovery_fit.rds> <query_pb_dir> <query_name> <out_prefix> [cts_file]
suppressPackageStartupMessages({
  library(cellprograms)
  library(rhdf5)
  library(clue)
})

args <- commandArgs(TRUE)
disc_path <- args[1]
pb_dir <- args[2]
qname <- args[3]
out_prefix <- args[4]
cts_filter <- if (length(args) >= 5 && nzchar(args[5])) readLines(args[5]) else NULL

.safe <- function(s) {
  z <- gsub("[^A-Za-z0-9_.-]+", "_", s)
  z <- sub("^_+", "", z)
  z <- sub("_+$", "", z)
  substr(z, 1, 80)
}

disc <- readRDS(disc_path)
stopifnot(!is.null(disc$canonicalization))

h5file <- file.path(pb_dir, "pseudobulk.h5")
genes <- h5read(h5file, "genes")
cts_all <- h5read(h5file, "cell_types")
cts <- if (is.null(cts_filter)) cts_all else intersect(cts_all, cts_filter)
meta <- read.csv(file.path(pb_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
if (!"observation_id" %in% names(meta)) meta$observation_id <- meta$donor_id

# keep the FULL query gene panel as columns: the per-ct projection intersects
# with each fit ct's own feature set, and pre-restricting to the union of
# discovery panels collapses per-ct coverage when panels differ across cts
shared_genes <- genes

# per-cell-type projection with tryCatch isolation
act_blocks <- list()
cov_rows <- list()
fail_rows <- list()
for (ct in cts) {
  g <- h5read(h5file, paste0("/", .safe(ct), "/counts"))
  d <- h5read(h5file, paste0("/", .safe(ct), "/donors"))
  if (nrow(g) == length(genes) && ncol(g) != length(genes)) g <- t(g)
  rownames(g) <- d
  colnames(g) <- genes
  lib <- rowSums(g)
  lib[lib == 0] <- 1
  m <- log1p(g[, shared_genes, drop = FALSE] / lib * 1e4)
  rownames(m) <- paste0(d, "|", ct)
  cell_types_ct <- stats::setNames(rep(ct, nrow(m)), rownames(m))
  donors_ct <- stats::setNames(d, rownames(m))
  out <- tryCatch(
    project_programs_to_cells(disc, m, cell_types_ct, donors = donors_ct,
                              sample_metadata = meta, min_coverage = 0.9),
    error = function(e) list(error = conditionMessage(e)))
  if (!is.null(out$error)) {
    fail_rows[[length(fail_rows) + 1]] <- data.frame(cell_type = ct, error = out$error)
    cat("FAIL:", ct, "-", out$error, "\n")
    next
  }
  act_blocks[[ct]] <- out$activity
  if (!is.null(out$coverage)) cov_rows[[length(cov_rows) + 1]] <- out$coverage
}

act <- NULL
all_ids <- unique(unlist(lapply(act_blocks, rownames)))
act <- matrix(NA_real_, nrow = length(all_ids), ncol = 0, dimnames = list(all_ids, NULL))
for (ct in names(act_blocks)) {
  blk <- act_blocks[[ct]]
  newcols <- setdiff(colnames(blk), colnames(act))
  if (length(newcols)) act <- cbind(act, matrix(NA_real_, nrow = length(all_ids), ncol = length(newcols), dimnames = list(all_ids, newcols)))
  act[rownames(blk), colnames(blk)] <- blk
}
if (length(cov_rows)) write.csv(do.call(rbind, cov_rows), paste0(out_prefix, "_coverage.csv"), row.names = FALSE)
if (length(fail_rows)) write.csv(do.call(rbind, fail_rows), paste0(out_prefix, "_failures.csv"), row.names = FALSE)

# concordance vs query refit per cell type
cons_rows <- list()
for (ct in intersect(names(act_blocks), disc$cell_types)) {
  cols <- grep(paste0("^", gsub("[][^$(){}+*?.|\\\\]", "\\\\\\0", ct), "_"), colnames(act), value = TRUE)
  if (!length(cols)) next
  ids_ct <- rownames(act_blocks[[ct]])
  sub <- act[ids_ct, cols, drop = FALSE]
  dvec <- sub("\\|[^|]*$", "", ids_ct)
  agg <- rowsum(sub, group = dvec, reorder = TRUE)
  refit_path <- file.path(pb_dir, paste0(tolower(qname), "_full_fit.rds"))
  if (file.exists(refit_path)) {
    rf <- readRDS(refit_path)
    if (!is.null(rf$scores[[ct]]) && ncol(rf$scores[[ct]]) > 0) {
      S <- rf$scores[[ct]]
      common <- intersect(rownames(agg), rownames(S))
      if (length(common) >= 10) {
        A <- agg[common, , drop = FALSE]
        B <- S[common, , drop = FALSE]
        ca <- cor(A, B)
        ca[is.na(ca)] <- 0
        if (min(dim(ca)) > 0) {
          if (nrow(ca) <= ncol(ca)) {
            mm <- clue::solve_LSAP(abs(ca), maximum = TRUE)
            pairs <- data.frame(proj = rownames(ca)[seq_along(mm)], refit = colnames(ca)[as.integer(mm)],
                                cor = ca[cbind(seq_along(mm), as.integer(mm))])
          } else {
            mm <- clue::solve_LSAP(t(abs(ca)), maximum = TRUE)
            pairs <- data.frame(proj = rownames(ca)[as.integer(mm)], refit = colnames(ca)[seq_along(mm)],
                                cor = ca[cbind(as.integer(mm), seq_along(mm))])
          }
          pairs$cell_type <- ct
          cons_rows[[length(cons_rows) + 1]] <- pairs
        }
      }
    }
  }
}

if (length(cons_rows)) {
  cons <- do.call(rbind, cons_rows)
  write.csv(cons, paste0(out_prefix, "_concordance.csv"), row.names = FALSE)
  cat("concordance written:", paste0(out_prefix, "_concordance.csv"),
      nrow(cons), "program pairs\n")
} else {
  cat("no concordance rows (refit fit not found or no shared cell types)\n")
}
saveRDS(act, paste0(out_prefix, "_activity.rds"))
cat("projection done for", qname, "\n")
