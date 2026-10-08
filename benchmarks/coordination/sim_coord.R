## sim_coord.R -- simulator for recurring programs across cell types.
## 10 cell types (ct01..ct10). Planted per dataset: R100 (all 10 cell types), R50 (5 random),
## R20 (2 random), 1-2 private programs per cell type, and an Institute nuisance on every cell type.
## mode crosses activity sharing with gene sharing:
##   act+genes : shared z, same genes (same loadings)      act+half : shared z, 50% of genes common
##   act_only  : shared z, cell-type-specific genes         genes_only: same genes, independent z per cell type
##   null      : no recurring programs (private + nuisance only)
## All planted gene sets are carved from one pool and are mutually disjoint.

COORD_CTS <- sprintf("ct%02d", 1:10)
COORD_BREADTH <- c(R100 = 10L, R50 = 5L, R20 = 2L)

gen_coord <- function(mode, amp, seed, act_cor = 1, G = getOption("ms.G", 4000L),
                      n_prog_genes = getOption("ms.npg", 50L), priv_amp = 1, n_donors = 100L) {
  stopifnot(mode %in% c("act+genes", "act+half", "act_only", "genes_only", "null"))
  set.seed(seed)
  cts <- COORD_CTS
  ids <- sprintf("d%03d", seq_len(n_donors))
  inst <- setNames(sample(c("I1", "I2", "I3"), n_donors, TRUE, prob = c(.6, .25, .15)), ids)
  bz <- setNames(as.numeric(scale(c(I1 = 0, I2 = 1, I3 = 2)[inst])), ids)
  ## donors per cell type: 4 large (all), 3 medium (75), 3 small (50; the small ones share a 30-donor core)
  sizes <- setNames(c(rep(100L, 4), rep(75L, 3), rep(50L, 3)), cts)[sample(10)]
  sizes <- sizes[cts]
  core <- sample(ids, 30)
  donors <- lapply(cts, function(ct) {
    n <- sizes[[ct]]
    if (n >= n_donors) ids
    else if (n == 50L) sort(c(core, sample(setdiff(ids, core), 20)))
    else sort(sample(ids, n))
  })
  names(donors) <- cts
  small <- cts[sizes == 50L]
  overlap <- outer(cts, cts, Vectorize(function(a, b) length(intersect(donors[[a]], donors[[b]]))))
  dimnames(overlap) <- list(cts, cts)

  genes <- paste0("g", seq_len(G))
  pool <- sample(G); ptr <- 0L
  take <- function(k) { stopifnot(ptr + k <= G); v <- pool[ptr + seq_len(k)]; ptr <<- ptr + k; v }
  bload <- rnorm(G) * (runif(G) < 0.5)

  ## planted recurring programs
  progs <- list()
  if (mode != "null") {
    for (nm in names(COORD_BREADTH)) {
      mem <- if (COORD_BREADTH[[nm]] == 10L) cts else sort(sample(cts, COORD_BREADTH[[nm]]))
      zs <- setNames(rnorm(n_donors), ids)
      gset <- w <- list()
      if (mode %in% c("act+genes", "genes_only")) {
        g <- take(n_prog_genes); wv <- rnorm(n_prog_genes, sd = amp)
        for (ct in mem) { gset[[ct]] <- g; w[[ct]] <- wv }
      } else if (mode == "act+half") {
        k <- n_prog_genes %/% 2L
        gc <- take(k); wc <- rnorm(k, sd = amp)
        for (ct in mem) { gs <- take(n_prog_genes - k); gset[[ct]] <- c(gc, gs); w[[ct]] <- c(wc, rnorm(length(gs), sd = amp)) }
      } else {   # act_only
        for (ct in mem) { gset[[ct]] <- take(n_prog_genes); w[[ct]] <- rnorm(n_prog_genes, sd = amp) }
      }
      z <- lapply(setNames(mem, mem), function(ct) {
        if (mode == "genes_only") setNames(rnorm(n_donors), ids)
        else if (act_cor >= 1) zs
        else act_cor * zs + sqrt(1 - act_cor^2) * setNames(rnorm(n_donors), ids)
      })
      progs[[nm]] <- list(id = nm, members = mem, z = z, genes = gset, w = w,
                          activity_true = mode != "genes_only", genes_true = mode != "act_only")
    }
  }
  priv <- lapply(setNames(cts, cts), function(ct) {
    lapply(seq_len(sample(1:2, 1)), function(k) list(g = take(n_prog_genes), w = rnorm(n_prog_genes, sd = priv_amp)))
  })

  pb <- list(); ncell <- list()
  for (ct in cts) {
    d <- donors[[ct]]; n <- length(d)
    nc <- setNames(exp(rnorm(n, log(500), 0.8)), d)
    sdi <- sqrt(stats::median(nc) / nc)
    Y <- matrix(rnorm(n * G), n, G, dimnames = list(d, genes)) * sdi
    Y <- Y + 1.5 * outer(bz[d], bload)
    for (p in priv[[ct]]) { zp <- rnorm(n); Y[, p$g] <- Y[, p$g] + outer(zp, p$w) }
    for (p in progs) if (ct %in% p$members) Y[, p$genes[[ct]]] <- Y[, p$genes[[ct]]] + outer(p$z[[ct]][d], p$w[[ct]])
    pb[[ct]] <- Y; ncell[[ct]] <- nc
  }
  meta <- data.frame(observation_id = ids, institute = factor(inst), stringsAsFactors = FALSE)
  list(x = as_cell_program_data(pb, sample_metadata = meta), ncell = ncell,
       truth = list(mode = mode, programs = progs, nuisance = bz, small = small, overlap = overlap),
       strata = setNames(as.character(inst), ids))
}
