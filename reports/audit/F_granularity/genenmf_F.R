# Workstream F (PILOT): GeneNMF arms on three input granularities.
# Package recipe (GeneNMF 0.9.2, Andreatta & Carmona): Seurat object,
# NormalizeData (log1p-CP10K) in the "data" slot, multiNMF = RcppML::nmf
# (MSE/Frobenius on log data, ONE run per K, no restarts) across K=2..15,
# then getMetaPrograms clusters top-gene sets across K into nMP meta-programs.
# hvg = the SAME top-2000 HVG list as the other arms (identical gene space).
# Pilot choices documented in TASKS.md: nMP=10 (package default), seed 20261017.
# All per-K w (genes x patterns) and h (patterns x obs) matrices are saved so
# the comparison script can build MP-level loadings/usages from memberships.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages(library(GeneNMF))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

root <- "/mnt/shared-workspace/shared/audit/F_granularity"
arms <- c("pseudobulk", "metacell", "singlecell")
SEED <- 20261017
KGRID <- 2:15
NMP <- 10

for (arm in arms) {
  d <- file.path(root, "gran", arm)
  X <- as.matrix(readMM(file.path(d, "counts.mtx")))
  meta <- read.delim(file.path(d, "obs_meta.tsv"), stringsAsFactors = FALSE)
  genes_raw <- readLines(file.path(d, "genes.tsv"))
  # Seurat v5 rewrites '_' to '-' in feature names, which breaks hvg subsetting
  # in getDataMatrix; sanitize BEFORE object creation and map back on output.
  genes <- gsub("_", "-", genes_raw)
  rownames(X) <- meta$id
  colnames(X) <- genes
  cat(sprintf("[%s] obs=%d genes=%d\n", arm, nrow(X), ncol(X)))

  obj <- CreateSeuratObject(counts = t(X),
                            meta.data = meta[, c("id", "donor"), drop = FALSE])
  obj <- NormalizeData(obj, verbose = FALSE)

  t0 <- Sys.time()
  nmf.res <- multiNMF(list(obj), k = KGRID, hvg = genes, seed = SEED)
  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # save all per-K solutions (restore original gene symbols in rownames)
  out_w <- NULL; out_h <- NULL
  for (nm in names(nmf.res)) {
    w <- nmf.res[[nm]]$w  # genes x patterns
    h <- nmf.res[[nm]]$h  # patterns x obs
    rownames(w) <- genes_raw
    colnames(w) <- paste0(nm, "_", colnames(w))
    rownames(h) <- paste0(nm, "_", rownames(h))
    out_w <- cbind(out_w, w)
    out_h <- rbind(out_h, h)  # patterns x obs: K rows per solution, K varies -> rbind
  }
  write.csv(out_w, sprintf("%s/F_genenmf_%s_allW.csv", root, arm))
  write.csv(out_h, sprintf("%s/F_genenmf_%s_allH.csv", root, arm))

  t1 <- Sys.time()
  mp <- tryCatch(
    getMetaPrograms(nmf.res, nMP = NMP),
    error = function(e) {
      cat(sprintf("[%s] getMetaPrograms failed (%s); using NaN-safe Jaccard fallback\n",
                  arm, conditionMessage(e)))
      pg <- getNMFgenes(nmf.res)
      pg <- pg[lengths(pg) > 0]
      if (length(pg) == 0) {
        return(list(markers.consensus = list(),
                    metrics = data.frame(metaProgram = character(0)),
                    composition = data.frame(metaProgram = character(0))))
      }
      allg <- unique(unlist(pg))
      bin <- t(vapply(pg, function(g) as.numeric(allg %in% g), numeric(length(allg))))
      inter <- bin %*% t(bin)
      sz <- rowSums(bin)
      J <- inter / pmax(outer(sz, sz, "+") - inter, 1e-12)
      diag(J) <- 1
      k_use <- min(NMP, nrow(bin))
      cl <- cutree(hclust(as.dist(1 - J), method = "average"), k = k_use)
      m <- length(unique(cl))
      markers <- lapply(seq_len(m), function(i) {
        members <- names(cl)[cl == i]
        cnt <- table(unlist(pg[members]))
        keep <- names(cnt)[cnt >= max(1, ceiling(length(members) / 2))]
        if (length(keep) == 0) keep <- names(cnt)[which.max(cnt)]
        keep
      })
      names(markers) <- paste0("MP", seq_len(m))
      list(markers.consensus = markers,
           metrics = data.frame(metaProgram = paste0("MP", seq_len(m)),
                                numberGenes = as.integer(lengths(markers)),
                                fallback = TRUE),
           composition = data.frame(metaProgram = paste0("MP", seq_len(m)),
                                    members = vapply(split(names(cl), cl),
                                                     paste, collapse = ";", "")))
    })
  t_sel <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("[%s] MPs=%d t_fit=%.0fs t_sel=%.0fs\n", arm,
              length(mp$markers.consensus), t_fit, t_sel))

  write.csv(data.frame(program = names(mp$markers.consensus),
                       genes = vapply(mp$markers.consensus, function(g) paste(g, collapse = ";"), ""),
                       stringsAsFactors = FALSE),
            sprintf("%s/F_genenmf_%s_mp_genes.csv", root, arm), row.names = FALSE)
  mp$metrics |> write.csv(sprintf("%s/F_genenmf_%s_mp_metrics.csv", root, arm))
  mp$composition |> write.csv(sprintf("%s/F_genenmf_%s_mp_composition.csv", root, arm), na = "")
  data.frame(arm = arm, n_mp = length(mp$markers.consensus), n_requested = NMP,
             k_grid = paste(KGRID, collapse = ","), seed = SEED,
             t_fit_s = t_fit, t_sel_s = t_sel) |>
    write.csv(sprintf("%s/F_genenmf_%s_summary.csv", root, arm), row.names = FALSE)
  rm(obj, nmf.res, out_w, out_h, mp); invisible(gc())
}
cat("GeneNMF F arms done\n")
