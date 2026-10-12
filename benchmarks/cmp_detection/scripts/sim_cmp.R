## sim_cmp.R -- frozen simulator for coordinated multicellular programs (CMPs).
## Two bases share one planting engine:
##   param : Gaussian pseudobulk (extends benchmarks/coordination/sim_coord.R), real gene symbols
##   semi  : real Stephenson pseudobulk, donors permuted within site independently per cell type
## Every planted latent is a "signal" with a type:
##   target  -> a CMP that should be called (shared donor activity, any genes)
##   decoy   -> must not be called (independent activity on same genes, ambient leak)
##   context -> target in marginal mode, decoy in within-context mode
## Scenarios: null, cmp_distinct, cmp_mixed, genes_only, context, ambient, lr, weak

CMP_RES <- "/mnt/shared-workspace/cmp/res"
CMP_SCEN <- c("null", "cmp_distinct", "cmp_mixed", "genes_only", "context", "ambient", "lr", "weak")
CMP_BREADTH <- c(R100 = 1, R50 = 0.5, R20 = 0.2)

.cmp_env <- new.env()
.cmp_res <- function() {
  if (is.null(.cmp_env$hall)) {
    h <- utils::read.csv(file.path(CMP_RES, "hallmark.csv"), stringsAsFactors = FALSE)
    .cmp_env$hall <- split(h$gene_symbol, h$gs_name)
    lr <- utils::read.csv(file.path(CMP_RES, "lr_consensus.csv"), stringsAsFactors = FALSE)
    .cmp_env$lr <- lr[!grepl("_", lr$ligand) & !grepl("_", lr$receptor), ]   # single-subunit pairs for planting
    .cmp_env$lr_all <- lr
  }
  .cmp_env
}

## ---- shared planting engine -----------------------------------------------------------
## state: list(cts, donors (per ct), pool (unused genes per ct or global), ...)
## take genes: n_theme from a Hallmark set + rest random, from the available pool
.take_genes <- function(avail, n, theme = NULL, frac_theme = 0.7) {
  g <- character()
  if (!is.null(theme)) {
    hs <- intersect(.cmp_res()$hall[[theme]], avail)
    g <- utils::head(sample(hs), round(frac_theme * n))
  }
  rest <- setdiff(avail, g)
  c(g, sample(rest, n - length(g)))
}

## spec of planted signals for a scenario. members drawn from cts.
.scenario_spec <- function(scenario, cts) {
  nct <- length(cts)
  br <- function(nm) if (nm == "R100") cts else sort(sample(cts, max(2L, round(CMP_BREADTH[[nm]] * nct))))
  mk <- function(id, members, genes_mode, type = "target", act_cor = 1, lr = FALSE)
    list(id = id, members = members, genes_mode = genes_mode, type = type, act_cor = act_cor, lr = lr)
  switch(scenario,
    null = list(),
    cmp_distinct = lapply(names(CMP_BREADTH), function(nm) mk(nm, br(nm), "distinct")),
    cmp_mixed = list(mk("R100", br("R100"), "same"), mk("R50", br("R50"), "half"), mk("R20", br("R20"), "distinct")),
    genes_only = lapply(names(CMP_BREADTH), function(nm) mk(nm, br(nm), "same", type = "decoy", act_cor = 0)),
    context = list(mk("CTX", sort(sample(cts, 4L)), "distinct", type = "context"), mk("R20", br("R20"), "distinct")),
    ambient = list(mk("R50", br("R50"), "distinct")),          # ambient decoy added separately
    lr = list(mk("R100", br("R100"), "distinct"), mk("R50", br("R50"), "distinct", lr = TRUE),
              mk("R20", br("R20"), "distinct", lr = TRUE)),
    weak = lapply(names(CMP_BREADTH), function(nm) mk(nm, br(nm), "distinct", act_cor = 0.7)),
    stop("unknown scenario ", scenario))
}

## plant signals into list of matrices Y (donor x gene), returns updated Y and truth
## wsd: function(ct) -> loading sd for that cell type; avail: function(ct) -> candidate genes
.plant <- function(Y, spec, ids, wsd, avail_fun, n_prog_genes, context_z = NULL) {
  cts <- names(Y)
  used <- setNames(lapply(cts, function(ct) character()), cts)
  themes <- names(.cmp_res()$hall)
  lrp <- .cmp_res()$lr
  sig <- list(); lr_truth <- list()
  for (s in spec) {
    zs <- setNames(rnorm(length(ids)), ids)
    z <- lapply(setNames(s$members, s$members), function(ct) {
      if (s$type == "context") { e <- setNames(rnorm(length(ids)), ids); sqrt(0.6) * context_z[ids] + sqrt(0.4) * e }
      else if (s$act_cor <= 0) setNames(rnorm(length(ids)), ids)
      else if (s$act_cor >= 1) zs
      else s$act_cor * zs + sqrt(1 - s$act_cor^2) * setNames(rnorm(length(ids)), ids)
    })
    genes <- list(); w <- list(); theme <- list()
    if (s$genes_mode %in% c("same", "half")) {
      common_avail <- Reduce(intersect, lapply(s$members, function(ct) setdiff(avail_fun(ct), used[[ct]])))
      th <- sample(themes, 1)
      k <- if (s$genes_mode == "same") n_prog_genes else n_prog_genes %/% 2L
      gc <- .take_genes(common_avail, k, th); wc <- rnorm(k)
    }
    for (ct in s$members) {
      av <- setdiff(avail_fun(ct), used[[ct]])
      if (s$genes_mode == "same") { g <- gc; wu <- wc; theme[[ct]] <- th }
      else if (s$genes_mode == "half") {
        th2 <- sample(themes, 1); gs <- .take_genes(setdiff(av, gc), n_prog_genes - k, th2)
        g <- c(gc, gs); wu <- c(wc, rnorm(length(gs))); theme[[ct]] <- th
      } else { th2 <- sample(themes, 1); g <- .take_genes(av, n_prog_genes, th2); wu <- rnorm(n_prog_genes); theme[[ct]] <- th2 }
      genes[[ct]] <- g; w[[ct]] <- wu
      used[[ct]] <- c(used[[ct]], g)
    }
    ## LR planting: ligand in sender's program, receptor in receiver(s)
    if (isTRUE(s$lr)) {
      snd <- s$members[1]; rcv <- s$members[-1][seq_len(min(2L, length(s$members) - 1L))]
      ok <- lrp[lrp$ligand %in% setdiff(avail_fun(snd), used[[snd]]), ]
      ok <- ok[vapply(ok$receptor, function(r) all(vapply(rcv, function(ct) r %in% setdiff(avail_fun(ct), used[[ct]]), TRUE)), TRUE), ]
      if (nrow(ok)) {
        p <- ok[sample(nrow(ok), 1), ]
        genes[[snd]][1] <- p$ligand; w[[snd]][1] <- 1.5 * sign(w[[snd]][1] + 1e-9)
        used[[snd]] <- c(used[[snd]], p$ligand)
        for (ct in rcv) { genes[[ct]][1] <- p$receptor; w[[ct]][1] <- 1.5 * sign(w[[ct]][1] + 1e-9); used[[ct]] <- c(used[[ct]], p$receptor) }
        lr_truth[[length(lr_truth) + 1L]] <- data.frame(signal = s$id, sender = snd, receiver = rcv,
                                                        ligand = p$ligand, receptor = p$receptor, stringsAsFactors = FALSE)
      }
    }
    for (ct in s$members) {
      d <- rownames(Y[[ct]])
      Y[[ct]][, genes[[ct]]] <- Y[[ct]][, genes[[ct]]] + outer(z[[ct]][d], w[[ct]] * wsd(ct))
    }
    sig[[s$id]] <- list(id = s$id, type = s$type, members = s$members, genes_mode = s$genes_mode,
                        z = lapply(setNames(names(z), names(z)), function(ct) z[[ct]][rownames(Y[[ct]])]),
                        genes = genes, w = w, theme = theme)
  }
  list(Y = Y, signals = sig, used = used,
       lr = if (length(lr_truth)) do.call(rbind, lr_truth) else
         data.frame(signal = character(), sender = character(), receiver = character(), ligand = character(), receptor = character()))
}

## ---- param base ---------------------------------------------------------------------------
.param_genes <- function(G = 4000L) {
  if (!is.null(.cmp_env$pgenes)) return(.cmp_env$pgenes)
  r <- .cmp_res()
  lrg <- unique(c(r$lr$ligand, r$lr$receptor))
  hg <- unique(unlist(r$hall))
  ## fixed gene universe drawn with a private seed; the caller's RNG stream is restored afterwards
  ## (bug fix: without restore, the first call in a fresh/forked process reset every dataset to seed 424242)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(424242)
  g <- unique(c(sample(lrg), sample(setdiff(hg, lrg))))[seq_len(G)]
  if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  .cmp_env$pgenes <- sort(g)
  .cmp_env$pgenes
}

gen_param <- function(scenario, seed, amp = 0.5, n_prog_genes = 50L, priv_amp = 1, n_donors = 100L) {
  set.seed(seed)
  genes <- .param_genes(); G <- length(genes)
  cts <- sprintf("ct%02d", 1:10)
  ids <- sprintf("d%03d", seq_len(n_donors))
  site <- setNames(sample(c("I1", "I2", "I3"), n_donors, TRUE, prob = c(.6, .25, .15)), ids)
  bz <- setNames(as.numeric(scale(c(I1 = 0, I2 = 1, I3 = 2)[site])), ids)
  sizes <- setNames(c(rep(100L, 4), rep(75L, 3), rep(50L, 3)), cts)[sample(10)][cts]
  core <- sample(ids, 30)
  donors <- lapply(setNames(cts, cts), function(ct) {
    n <- sizes[[ct]]
    if (n >= n_donors) ids else if (n == 50L) sort(c(core, sample(setdiff(ids, core), 20))) else sort(sample(ids, n))
  })
  bload <- setNames(rnorm(G) * (runif(G) < 0.5), genes)
  Y <- list(); ncell <- list()
  for (ct in cts) {
    d <- donors[[ct]]; n <- length(d)
    nc <- setNames(exp(rnorm(n, log(500), 0.8)), d)
    Y[[ct]] <- matrix(rnorm(n * G), n, G, dimnames = list(d, genes)) * sqrt(stats::median(nc) / nc) +
      1.5 * outer(bz[d], bload)
    ncell[[ct]] <- nc
  }
  ctx <- if (scenario == "context") { c0 <- setNames(as.numeric(runif(n_donors) < 0.4), ids); setNames(as.numeric(scale(c0)), ids) } else NULL
  spec <- .scenario_spec(scenario, cts)
  pl <- .plant(Y, spec, ids, wsd = function(ct) amp, avail_fun = function(ct) genes, n_prog_genes, context_z = ctx)
  Y <- pl$Y; used <- pl$used; signals <- pl$signals
  ## private programs (themed), and decoy LR genes placed in private programs of unrelated cell types
  themes <- names(.cmp_res()$hall)
  priv <- list()
  for (ct in cts) for (k in seq_len(sample(1:2, 1))) {
    g <- .take_genes(setdiff(genes, used[[ct]]), n_prog_genes, sample(themes, 1)); used[[ct]] <- c(used[[ct]], g)
    zp <- rnorm(nrow(Y[[ct]])); wv <- rnorm(n_prog_genes, sd = priv_amp)
    Y[[ct]][, g] <- Y[[ct]][, g] + outer(zp, wv)
    priv[[paste(ct, k)]] <- list(ct = ct, genes = g, z = setNames(zp, rownames(Y[[ct]])))
  }
  ## ambient decoy: sender marker profile leaks into every other cell type with a donor-level fraction
  if (scenario == "ambient") {
    snd <- names(sizes)[which.max(sizes)]
    mk <- sample(setdiff(genes, unlist(lapply(signals, function(s) unlist(s$genes)))), 150)
    m <- abs(rnorm(150)) * 2
    f <- setNames(as.numeric(scale(exp(rnorm(n_donors, 0, 0.8)))), ids)
    rcv <- setdiff(cts, snd)
    for (ct in rcv) { d <- rownames(Y[[ct]]); Y[[ct]][, mk] <- Y[[ct]][, mk] + 0.6 * amp * outer(f[d], m) }
    signals$AMB <- list(id = "AMB", type = "decoy", members = rcv, genes_mode = "same",
                        z = lapply(setNames(rcv, rcv), function(ct) f[rownames(Y[[ct]])]),
                        genes = setNames(rep(list(mk), length(rcv)), rcv), sender = snd)
  }
  meta <- data.frame(observation_id = ids, institute = factor(site), stringsAsFactors = FALSE)
  if (!is.null(ctx)) meta$context <- ctx
  list(x = as_cell_program_data(Y, sample_metadata = meta), ncell = ncell,
       strata = setNames(as.character(site), ids), context = ctx, site_z = bz,
       truth = list(base = "param", scenario = scenario, seed = seed, signals = signals, lr = pl$lr,
                    private = priv, sizes = sizes))
}

## ---- semi base (Stephenson) -------------------------------------------------------------
SEMI_DIR <- "/mnt/shared-workspace/cmp/data/stephenson_pb"
.semi_load <- function() {
  if (!is.null(.cmp_env$semi)) return(.cmp_env$semi)
  .cmp_env$semi <- readRDS(file.path(SEMI_DIR, "semi_base.rds"))   # list(ct = list(Y logCPM, site, ncell, msd, cpm_mean), ...)
  .cmp_env$semi
}

gen_semi <- function(scenario, seed, amp = 0.5, n_prog_genes = 50L) {
  b <- .semi_load()
  set.seed(seed)
  cts <- names(b$cts)
  Y <- list(); ncell <- list(); site_of <- c()
  for (ct in cts) {
    M <- b$cts[[ct]]$Y; lab <- rownames(M); st <- b$cts[[ct]]$site; new <- lab
    for (i in unique(st)) { idx <- which(st == i); new[idx] <- lab[idx][sample.int(length(idx))] }
    rownames(M) <- new; nc <- b$cts[[ct]]$ncell; names(nc) <- new
    Y[[ct]] <- M; ncell[[ct]] <- nc; site_of[new] <- as.character(st)
  }
  ids <- sort(unique(unlist(lapply(Y, rownames))))
  ctx <- if (scenario == "context") { c0 <- setNames(as.numeric(runif(length(ids)) < 0.4), ids); setNames(as.numeric(scale(c0)), ids) } else NULL
  spec <- .scenario_spec(scenario, cts)
  pl <- .plant(Y, spec, ids, wsd = function(ct) amp * b$cts[[ct]]$msd,
               avail_fun = function(ct) colnames(b$cts[[ct]]$Y), n_prog_genes, context_z = ctx)
  Y <- pl$Y; signals <- pl$signals
  if (scenario == "ambient") {
    snd <- b$ambient_sender
    prof <- b$cts[[snd]]$cpm_mean
    f <- setNames(0.03 * exp(rnorm(length(ids), 0, 0.8)), ids)
    rcv <- setdiff(cts, snd)
    for (ct in rcv) {
      d <- rownames(Y[[ct]]); g <- intersect(colnames(Y[[ct]]), names(prof))
      cpm <- 2^Y[[ct]][, g] - 1
      Y[[ct]][, g] <- log2((1 - f[d]) * cpm + outer(f[d], prof[g]) + 1)
    }
    top <- names(sort(prof, decreasing = TRUE))[1:150]
    signals$AMB <- list(id = "AMB", type = "decoy", members = rcv, genes_mode = "same",
                        z = lapply(setNames(rcv, rcv), function(ct) as.numeric(scale(f[rownames(Y[[ct]])])) |> setNames(rownames(Y[[ct]]))),
                        genes = setNames(rep(list(top), length(rcv)), rcv), sender = snd)
  }
  meta <- data.frame(observation_id = ids, institute = factor(site_of[ids]), stringsAsFactors = FALSE)
  if (!is.null(ctx)) meta$context <- ctx
  list(x = as_cell_program_data(Y, sample_metadata = meta), ncell = ncell,
       strata = setNames(as.character(site_of[ids]), ids), context = ctx,
       truth = list(base = "semi", scenario = scenario, seed = seed, signals = signals, lr = pl$lr))
}

gen_cmp <- function(base, scenario, seed, amp) {
  off <- (match(base, c("param", "semi")) - 1L) * 5e6 + match(scenario, CMP_SCEN) * 1e5
  if (base == "param") gen_param(scenario, off + seed, amp) else gen_semi(scenario, off + seed, amp)
}
