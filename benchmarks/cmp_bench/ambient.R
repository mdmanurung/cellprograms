## ambient.R -- F8: ambient / cross-cell-type leakage score per node.
## Inputs: per-cell-type mean expression profiles (CPM-like, named numeric vectors). For node (ct, k) the
## ambient score is max over other cell types s of cor(loading_k, log2(prof_s+1) - log2(prof_ct+1)) over the
## node's genes (positive direction only, after sign-aligning the loading so that the score is maximal).
## Real data: profiles = per-cell-type pseudobulk mean CPM. Benchmark: semi uses the real Stephenson means;
## param uses synthetic profiles (lognormal baseline; the largest cell type carries 150 marker genes, which in
## the ambient scenario are exactly the leaked genes).
SEMI_BASE_CACHE <- new.env()
get_profiles <- function(o) {
  if (o$truth$base == "semi") {
    if (is.null(SEMI_BASE_CACHE$p)) {
      b <- readRDS(file.path("/mnt/shared-workspace/cmp/data/stephenson_pb/semi_base.rds"))
      SEMI_BASE_CACHE$p <- lapply(b$cts, `[[`, "cpm_mean")
    }
    return(SEMI_BASE_CACHE$p)
  }
  genes <- .param_genes(); cts <- names(o$truth$sizes)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(o$truth$seed + 777L)
  base <- exp(stats::rnorm(length(genes), 2, 1))
  snd <- names(o$truth$sizes)[which.max(o$truth$sizes)]
  mk <- if (!is.null(o$truth$signals$AMB)) o$truth$signals$AMB$genes[[1]] else sample(genes, 150)
  prof <- lapply(stats::setNames(cts, cts), function(ct) stats::setNames(base * exp(stats::rnorm(length(genes), 0, 0.5)), genes))
  prof[[snd]][mk] <- prof[[snd]][mk] * 2^(2 + abs(stats::rnorm(length(mk))))
  if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  prof
}
ambient_scores <- function(fit, profiles) {
  lp <- lapply(profiles, function(p) log2(p + 1))
  do.call(rbind, lapply(fit$cell_types, function(ct) {
    W <- fit$loadings[[ct]]; if (is.null(W) || !ncol(W) || is.null(lp[[ct]])) return(NULL)
    best <- rep(0, ncol(W)); src <- rep(NA_character_, ncol(W))
    for (s in setdiff(names(lp), ct)) {
      g <- intersect(rownames(W), intersect(names(lp[[s]]), names(lp[[ct]])))
      if (length(g) < 50) next
      spec <- lp[[s]][g] - lp[[ct]][g]
      r <- abs(suppressWarnings(stats::cor(W[g, , drop = FALSE], spec))); r[!is.finite(r)] <- 0
      up <- r > best; best[up] <- r[up]; src[up] <- s
    }
    data.frame(cell_type = ct, program_id = colnames(W), amb = best, source = src, stringsAsFactors = FALSE)
  }))
}
