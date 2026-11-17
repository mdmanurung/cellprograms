## sim_coord.R -- simulator for recurring programs across cell types.
## 10 cell types (ct01..ct10). Modes:
##   act+genes : shared z, same genes (same loadings)      act+half : shared z, 50% of genes common
##   act_only  : shared z, cell-type-specific genes         genes_only: same genes, independent z per cell type
##   priv_corr : NO recurring programs; 3 pairs of cell types get private programs with
##               correlated donor activity (cor_priv) but disjoint genes (the S5 trap:
##               activity-correlated but genuinely distinct programs must NOT be called shared)
##   null      : no recurring programs (private + nuisance only)
## breadth: "all" co-plants R100+R50+R20; "R100"/"R50"/"R20" plant that program alone;
##   "R20small" plants an R20 whose members are forced into two small cell types.
## All planted gene sets are carved from one pool and are mutually disjoint.
## Trap modes (added 2026-11-17; no recurring programs, so any multi-cell-type call is spurious).
##   donor_deletion: a set of n_del_genes genes is absent (shifted down by del_depth) in a fraction
##                   del_frac of donors, in every cell type -- a genotype-like binary donor axis.
##   ig_tcr        : in n_clone_cts cell types, a set of n_clone_genes genes is large in a fraction
##                   clone_frac of that cell type's donors (drawn independently per cell type),
##                   with gene weights of sd clone_mult * amp -- clonal Ig/TCR-like expression.
## These arguments have no defaults: their values are impactful design choices (prereg-coord-v2).

## `.residualize()` was an internal helper of R/ebmf.R until 0.3.0 and is gone from the package; this
## benchmark code still calls it. Same definition as before, on top of the package's .cov_design().
if (!exists(".residualize", mode = "function")) {
  .residualize <- function(Y, meta, covariates, ct) {
    X <- .cov_design(rownames(Y), meta, covariates, ct)
    if (is.null(X)) return(Y)
    R <- stats::lm.fit(X, Y)$residuals
    dimnames(R) <- dimnames(Y)
    R
  }
}

COORD_CTS <- sprintf("ct%02d", 1:10)
COORD_BREADTH <- c(R100 = 10L, R50 = 5L, R20 = 2L)

gen_coord <- function(mode, amp, seed, act_cor = 1, G = getOption("ms.G", 4000L),
                      n_prog_genes = getOption("ms.npg", 50L), priv_amp = 1, n_donors = 100L,
                      breadth = c("all", "R100", "R50", "R20", "R20small"),
                      cor_priv = 0.7, n_decoy_pairs = 3L, private_per_ct = c(1L, 2L),
                      n_priv_genes = n_prog_genes,
                      del_frac = NULL, n_del_genes = NULL, del_depth = NULL,
                      clone_frac = NULL, clone_mult = NULL, n_clone_cts = NULL, n_clone_genes = NULL) {
  stopifnot(mode %in% c("act+genes", "act+half", "act_only", "genes_only", "priv_corr", "null",
                        "donor_deletion", "ig_tcr"))
  if (mode == "donor_deletion") {
    stopifnot(!is.null(del_frac), !is.null(n_del_genes), !is.null(del_depth))
  }
  if (mode == "ig_tcr") {
    stopifnot(!is.null(clone_frac), !is.null(clone_mult), !is.null(n_clone_cts), !is.null(n_clone_genes))
  }
  breadth <- match.arg(breadth)
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

  ## which recurring programs to plant
  plant <- switch(breadth,
                  all  = names(COORD_BREADTH),
                  R100 = "R100", R50 = "R50", R20 = "R20", R20small = "R20")
  if (mode %in% c("null", "priv_corr", "donor_deletion", "ig_tcr")) plant <- character(0)

  ## planted recurring programs
  progs <- list()
  for (nm in plant) {
    mem <- if (breadth == "R20small") {
      sort(sample(small, 2L))
    } else if (COORD_BREADTH[[nm]] == 10L) {
      cts
    } else {
      sort(sample(cts, COORD_BREADTH[[nm]]))
    }
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

  ## correlated private pairs (priv_corr): distinct genes, correlated activity.
  ## Truth: these are DISTINCT programs; a multi-cell-type call over a pair is a trap failure.
  decoy_pairs <- list()
  if (mode == "priv_corr") {
    pair_cts <- sample(cts, 2L * n_decoy_pairs)
    for (j in seq_len(n_decoy_pairs)) {
      a <- pair_cts[2L * j - 1L]; b <- pair_cts[2L * j]
      zb <- setNames(rnorm(n_donors), ids)
      za <- zb
      zb2 <- cor_priv * zb + sqrt(1 - cor_priv^2) * setNames(rnorm(n_donors), ids)
      decoy_pairs[[j]] <- list(
        id = sprintf("DP%d", j), members = c(a, b),
        z = setNames(list(za, zb2), c(a, b)),
        genes = setNames(list(take(n_prog_genes), take(n_prog_genes)), c(a, b)),
        w = setNames(list(rnorm(n_prog_genes, sd = amp), rnorm(n_prog_genes, sd = amp)), c(a, b))
      )
    }
  }

  ## private_per_ct / n_priv_genes (post-prereg high-K scenario): number of
  ## private programs per cell type and their gene-set size. Defaults
  ## (c(1L, 2L) drawn per ct; n_priv_genes = n_prog_genes) reproduce the
  ## frozen main-run simulator exactly.
  priv <- lapply(setNames(cts, cts), function(ct) {
    lapply(seq_len(sample(private_per_ct, 1)),
           function(k) list(g = take(n_priv_genes), w = rnorm(n_priv_genes, sd = priv_amp)))
  })

  ## trap structures (new modes only; their draws come after every draw the older modes make)
  traps <- list()
  if (mode == "donor_deletion") {
    traps$donor_deletion <- list(genes = take(n_del_genes),
                                 donors = sort(sample(ids, round(del_frac * n_donors))),
                                 depth = del_depth)
  }
  if (mode == "ig_tcr") {
    chosen <- sort(sample(cts, n_clone_cts))
    traps$ig_tcr <- lapply(setNames(chosen, chosen), function(ct) {
      list(genes = take(n_clone_genes),
           donors = sort(sample(donors[[ct]], max(1L, round(clone_frac * length(donors[[ct]]))))),
           w = rnorm(n_clone_genes, sd = clone_mult * amp))
    })
  }

  pb <- list(); ncell <- list()
  priv_truth <- setNames(vector("list", length(cts)), cts)
  for (ct in cts) {
    d <- donors[[ct]]; n <- length(d)
    nc <- setNames(exp(rnorm(n, log(500), 0.8)), d)
    sdi <- sqrt(stats::median(nc) / nc)
    Y <- matrix(rnorm(n * G), n, G, dimnames = list(d, genes)) * sdi
    Y <- Y + 1.5 * outer(bz[d], bload)
    pt <- list()
    for (p in priv[[ct]]) {
      zp <- setNames(rnorm(n), d)
      Y[, p$g] <- Y[, p$g] + outer(zp, p$w)
      pt[[length(pt) + 1L]] <- list(z = zp, genes = p$g, w = p$w)
    }
    for (p in progs) if (ct %in% p$members) Y[, p$genes[[ct]]] <- Y[, p$genes[[ct]]] + outer(p$z[[ct]][d], p$w[[ct]])
    for (dp in decoy_pairs) if (ct %in% dp$members) {
      Y[, dp$genes[[ct]]] <- Y[, dp$genes[[ct]]] + outer(dp$z[[ct]][d], dp$w[[ct]])
      pt[[length(pt) + 1L]] <- list(z = dp$z[[ct]][d], genes = dp$genes[[ct]], w = dp$w[[ct]],
                                    decoy_pair = dp$id)
    }
    if (!is.null(traps$donor_deletion)) {          # genes absent in the deleted donors, every cell type
      tp <- traps$donor_deletion; rows <- which(d %in% tp$donors)
      if (length(rows)) Y[rows, tp$genes] <- Y[rows, tp$genes] - tp$depth
    }
    if (!is.null(traps$ig_tcr[[ct]])) {            # clonal donors carry the gene set at high weight
      ti <- traps$ig_tcr[[ct]]
      Y[, ti$genes] <- Y[, ti$genes] + outer(as.numeric(d %in% ti$donors), ti$w)
    }
    priv_truth[[ct]] <- pt
    pb[[ct]] <- Y; ncell[[ct]] <- nc
  }
  meta <- data.frame(observation_id = ids, institute = factor(inst), stringsAsFactors = FALSE)
  out <- list(x = as_cell_program_data(pb, sample_metadata = meta), ncell = ncell,
              truth = list(mode = mode, breadth = breadth, programs = progs, decoy_pairs = decoy_pairs,
                           privates = priv_truth, nuisance = bz, small = small, overlap = overlap,
                           donors = donors),
              strata = setNames(as.character(inst), ids))
  if (length(traps)) out$truth$traps <- traps    # only the trap modes carry this element
  out
}

## Adjusted expression matrices without stage-1 fits: residualize on the
## measured nuisance covariate and gene-center, mirroring the preregistered
## pre-stage-1 adjustment position. Used by eb_cp/block_E inputs and the oracle.
adjusted_matrices <- function(x, covariates = "institute") {
  out <- list()
  for (ct in x$cell_types) {
    Y <- as.matrix(x$pseudobulk[[ct]])
    if (!is.null(covariates)) Y <- .residualize(Y, x$sample_metadata, covariates, ct)
    Y <- scale(Y, center = TRUE, scale = FALSE)
    cs <- colSums(abs(Y))
    Y <- Y[, cs > 1e-10 * max(cs, 1), drop = FALSE]
    out[[ct]] <- Y
  }
  out
}
