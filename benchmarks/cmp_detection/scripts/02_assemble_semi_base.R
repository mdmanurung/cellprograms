#!/usr/bin/env Rscript
## 02_assemble_semi_base.R -- assemble the rebuilt Stephenson substrate into
## semi_base.rds (the contract real_data.R consumes):
##   list(cts = list(<type> = list(Y = donor x gene counts, cpm_mean = gene
##        vector)), donor_meta = data.frame(donor_id, site, sex,
##        Status_on_day_collection_summary))
## Usage: Rscript 02_assemble_semi_base.R [semi_dir]

semi_dir <- commandArgs(TRUE)[1] %||% file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "data", "stephenson_pb")
`%||%` <- function(a, b) if (is.null(a)) b else a
semi_dir <- commandArgs(TRUE)[1] %||% semi_dir
suppressPackageStartupMessages(library(Matrix))

types <- sort(list.dirs(semi_dir, recursive = FALSE, full.names = FALSE))
stopifnot(length(types) >= 2L)
donor_meta <- read.csv(file.path(semi_dir, "donor_meta.csv"), stringsAsFactors = FALSE)
cts <- lapply(types, function(ty) {
  td <- file.path(semi_dir, ty)
  Y <- readMM(gzfile(file.path(td, "Y.mtx.gz")))
  Y <- as.matrix(Y)
  dimnames(Y) <- list(readLines(file.path(td, "donors.txt")),
                      as.character(read.csv(file.path(td, "cpm_mean.csv"))$gene))
  ## Drop CITE-seq ADT features (AB_*): they are protein-derived (centered,
  ## negative-valued), not RNA counts. The lost original semi_base excluded them
  ## (RNA-only gene universe); the rebuild restores that convention. Verified:
  ## all negative values in the rebuilt Y sit in AB_* columns.
  ab <- grep("^AB_", colnames(Y))
  if (length(ab)) Y <- Y[, -ab, drop = FALSE]
  list(Y = Y, cpm_mean = setNames(read.csv(file.path(td, "cpm_mean.csv"))$cpm_mean,
                                  colnames(Y)),
       msd = stats::median(apply(Y, 2, stats::sd), na.rm = TRUE),
       site = local({ dm <- donor_meta; s <- dm$site[match(rownames(Y), dm$donor_id)];
         setNames(as.character(s), rownames(Y)) }),
       ncell = local({ ncf <- file.path(td, "ncell.csv");
         if (file.exists(ncf)) { ncv <- read.csv(ncf); v <- setNames(as.integer(ncv$ncell), as.character(ncv$donor));
           v[rownames(Y)] %||% setNames(rep(1L, nrow(Y)), rownames(Y)) } else setNames(rep(1L, nrow(Y)), rownames(Y)) }))
})
names(cts) <- types
stopifnot(all(c("donor_id", "site", "sex", "Status_on_day_collection_summary") %in% names(donor_meta)))

## ambient_sender: the largest cell type (donor count), matching the param-base
## convention (sim_cmp.R param branch: names(sizes)[which.max(sizes)]). The lost
## original builder set this field; the rebuild restores it deterministically.
sizes <- vapply(cts, function(x) nrow(x$Y), integer(1))
ambient_sender <- names(sizes)[which.max(sizes)]
out <- list(cts = cts, donor_meta = donor_meta, ambient_sender = ambient_sender)
f <- file.path(semi_dir, "semi_base.rds")
tmp <- tempfile(fileext = ".rds"); saveRDS(out, tmp)
ok <- system2("cp", c(shQuote(tmp), shQuote(f))) == 0L; unlink(tmp)
stopifnot(ok && file.size(f) > 0)
cat("semi_base.rds written:", f, "\n")
for (ty in types) cat(sprintf("  %s: %d donors x %d genes\n", ty,
                              nrow(cts[[ty]]$Y), ncol(cts[[ty]]$Y)))
