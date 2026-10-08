## sim_param.R -- parametric generator for the model-selection benchmark.
## 4 cell types: A, B (100 donors), C (d001-d050), D (d021-d070; 30 donors shared with C).
## Batch (Institute) shifts the same genes in every cell type; n_cells is lognormal and
## noise sd_i = sqrt(median(n)/n_i) (exactly the w_S model; see PREREG).
## scenario: S0 null | S1 shared A-B | S2 shared A-C | S3 shared C-D | S7 same genes A-B, independent z.

.SCEN_PAIR <- list(S1 = c("A", "B"), S2 = c("A", "C"), S3 = c("C", "D"), S7 = c("A", "B"))

gen_param <- function(scenario, amp, seed, G = getOption("ms.G", 4000L), n_prog_genes = 50L,
                      n_priv = 2L, priv_amp = 1) {
  set.seed(seed)
  ids <- sprintf("d%03d", 1:100)
  inst <- sample(c("I1", "I2", "I3"), 100, TRUE, prob = c(.6, .25, .15))
  names(inst) <- ids
  bz <- setNames(as.numeric(scale(c(I1 = 0, I2 = 1, I3 = 2)[inst])), ids)
  donors <- list(A = ids, B = ids, C = ids[1:50], D = ids[21:70])
  genes <- paste0("g", seq_len(G))
  bload <- rnorm(G) * (runif(G) < 0.5)            # shared by all cell types
  pair <- .SCEN_PAIR[[scenario]]
  zs <- setNames(rnorm(100), ids)                  # shared activity over all donors
  w_unit <- rnorm(n_prog_genes); g_same <- sample(G, n_prog_genes)   # for S7
  pb <- list(); ncell <- list(); z_truth <- list()
  for (ct in names(donors)) {
    d <- donors[[ct]]; n <- length(d)
    nc <- setNames(exp(rnorm(n, log(500), 0.8)), d)
    sdi <- sqrt(median(nc) / nc)
    Y <- matrix(rnorm(n * G), n, G, dimnames = list(d, genes)) * sdi
    Y <- Y + 1.5 * outer(bz[d], bload)
    for (k in seq_len(n_priv)) {
      w <- numeric(G); w[sample(G, n_prog_genes)] <- rnorm(n_prog_genes, sd = priv_amp)
      Y <- Y + outer(rnorm(n), w)
    }
    if (!is.null(pair) && ct %in% pair) {
      z <- if (scenario == "S7") setNames(rnorm(100), ids) else zs
      w <- numeric(G)
      if (scenario == "S7") w[g_same] <- w_unit * amp
      else w[sample(G, n_prog_genes)] <- rnorm(n_prog_genes, sd = amp)
      Y <- Y + outer(z[d], w)
      z_truth[[ct]] <- z[d]
    }
    pb[[ct]] <- Y; ncell[[ct]] <- nc
  }
  meta <- data.frame(observation_id = ids, institute = factor(inst), stringsAsFactors = FALSE)
  list(x = as_cell_program_data(pb, sample_metadata = meta), ncell = ncell,
       truth = list(scenario = scenario, pair = pair, z = z_truth),
       strata = setNames(as.character(inst), ids))
}
