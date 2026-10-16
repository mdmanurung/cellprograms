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
                      n_prog_genes = getOption("ms.npg", 50L), priv_amp = 1, n_donors = 100L,
                      artifacts = character(0)) {
  stopifnot(mode %in% c("act+genes", "act+half", "act_only", "genes_only", "null"))
  stopifnot(all(artifacts %in% c("vdj", "deletion")))
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

  ## Artifact-shaped decoys (HANDOFF 2026-10-11 adoption items). Neither is a target: clusters
  ## aligned with them are artifact calls, reported like decoys, never in the F1.
  ## - vdj: an Ig/TCR-like program -- a donor-composition latent expressed on a small receptor-gene
  ##   pool (renamed to IGLV*/TRBV* so the package's artifact_genes() regex flags it) in a random
  ##   subset of 3-5 cell types, same genes across members (the real shape).
  ## - deletion: 1-2 donors shifted on a shared gene set in EVERY cell type (a donor property);
  ##   the "activity" is the donor indicator, the exact shape of a low-quality-donor artifact.
  art <- list()
  vdj_genes <- NULL
  if ("vdj" %in% artifacts) {
    vdj_genes <- take(40L)
    mem_v <- sort(sample(cts, sample(3:5, 1)))
    zv <- setNames(rnorm(n_donors), ids)
    wv <- rnorm(40L, sd = amp)
    art[["VDJ"]] <- list(id = "VDJ", type = "artifact_vdj", members = mem_v,
                         z = setNames(lapply(mem_v, function(ct) zv), mem_v),
                         genes = setNames(rep(list(vdj_genes), length(mem_v)), mem_v),
                         w = setNames(rep(list(wv), length(mem_v)), mem_v),
                         activity_true = TRUE, genes_true = TRUE)
  }
  if ("deletion" %in% artifacts) {
    d_del <- sample(ids, sample(1:2, 1))
    g_del <- take(100L)
    w_del <- rnorm(100L, sd = amp)
    zd <- setNames(as.numeric(ids %in% d_del), ids)
    art[["DEL"]] <- list(id = "DEL", type = "artifact_deletion", members = cts,
                         z = setNames(lapply(cts, function(ct) zd), cts),
                         genes = setNames(rep(list(g_del), length(cts)), cts),
                         w = setNames(rep(list(w_del), length(cts)), cts),
                         activity_true = TRUE, genes_true = TRUE)
  }

  pb <- list(); ncell <- list()
  for (ct in cts) {
    d <- donors[[ct]]; n <- length(d)
    nc <- setNames(exp(rnorm(n, log(500), 0.8)), d)
    sdi <- sqrt(stats::median(nc) / nc)
    Y <- matrix(rnorm(n * G), n, G, dimnames = list(d, genes)) * sdi
    Y <- Y + 1.5 * outer(bz[d], bload)
    for (p in priv[[ct]]) { zp <- rnorm(n); Y[, p$g] <- Y[, p$g] + outer(zp, p$w) }
    for (p in progs) if (ct %in% p$members) Y[, p$genes[[ct]]] <- Y[, p$genes[[ct]]] + outer(p$z[[ct]][d], p$w[[ct]])
    for (p in art) if (ct %in% p$members) Y[, p$genes[[ct]]] <- Y[, p$genes[[ct]]] + outer(p$z[[ct]][d], p$w[[ct]])
    pb[[ct]] <- Y; ncell[[ct]] <- nc
  }
  if (!is.null(vdj_genes)) {
    vdj_names <- c(paste0("IGLV", seq_len(20L)), paste0("TRBV", seq_len(20L)))
    for (ct in cts) colnames(pb[[ct]])[vdj_genes] <- vdj_names
  }
  meta <- data.frame(observation_id = ids, institute = factor(inst), stringsAsFactors = FALSE)
  list(x = as_cell_program_data(pb, sample_metadata = meta), ncell = ncell,
       truth = list(mode = mode, programs = c(progs, art), nuisance = bz, small = small, overlap = overlap),
       strata = setNames(as.character(inst), ids))
}

## Count-level variant (HANDOFF: "count-level simulator (current is Gaussian)"). Same latent
## structure and truth as gen_coord, but expression is realized as Poisson counts per donor
## (library size ncell, relative expression exp(Y) normalized per donor) and the fit input is
## log2 CPM. Stated limit: Poisson aggregation only -- no NB overdispersion arm.
gen_coord_counts <- function(mode, amp, seed, ...) {
  d <- gen_coord(mode, amp, seed, ...)
  set.seed(seed + 1L)
  pb_c <- stats::setNames(lapply(names(d$x$pseudobulk), function(ct) {
    Y <- d$x$pseudobulk[[ct]]
    nc <- d$ncell[[ct]][rownames(Y)]
    rho <- exp(Y)
    lambda <- rho / rowSums(rho) * nc
    matrix(stats::rpois(length(lambda), lambda), nrow(Y), ncol(Y), dimnames = dimnames(Y))
  }), names(d$x$pseudobulk))
  pb_lcpm <- lapply(pb_c, function(C) log2(sweep(C, 1L, rowSums(C), "/") * 1e6 + 1))
  d$x <- as_cell_program_data(pb_lcpm, sample_metadata = d$x$sample_metadata)
  d$counts <- pb_c
  d
}
