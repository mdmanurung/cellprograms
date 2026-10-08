## sim_semi.R -- semi-synthetic generator on real COMBAT pseudobulk.
## A=CD4, B=cMono, C=GDT, D=DP. Per rep, donors are permuted within Institute independently
## per cell type (rows keep their n_cells), which removes real cross-cell-type biology but
## keeps Institute structure, real noise and real heteroscedasticity. Programs are planted
## on 50 genes drawn from the full matrix: w ~ N(0, (amp * median gene SD)^2), z ~ N(0,1).

.semi_env <- new.env()
.semi_load <- function(root) {
  if (!is.null(.semi_env$d)) return(.semi_env$d)
  dir <- file.path(root, "benchmarks", "biomni_replication", "data", "combat")
  meta <- utils::read.csv(file.path(dir, "donor_meta.csv"), stringsAsFactors = FALSE)
  pbc <- utils::read.csv(file.path(dir, "pb_coldata.csv"), stringsAsFactors = FALSE)
  cts <- c(A = "CD4", B = "cMono", C = "GDT", D = "DP")
  d <- lapply(cts, function(ct) {
    Y <- as.matrix(utils::read.csv(file.path(dir, "mat", paste0("ct_", ct, ".csv")),
                                   row.names = 1, check.names = FALSE))
    inst <- meta$Institute[match(rownames(Y), meta$donor_id)]
    ok <- !is.na(inst)
    Y <- Y[ok, , drop = FALSE]; inst <- inst[ok]
    p <- pbc[pbc$cell_type == ct, ]
    list(Y = Y, inst = setNames(inst, rownames(Y)),
         ncell = setNames(p$n_cells[match(rownames(Y), p$donor_id)], rownames(Y)),
         msd = stats::median(apply(Y, 2, stats::sd)))
  })
  .semi_env$d <- d
  d
}

gen_semi <- function(scenario, amp, seed, root, n_prog_genes = 50L) {
  d <- .semi_load(root)
  set.seed(seed)
  pair <- .SCEN_PAIR[[scenario]]
  pb <- list(); ncell <- list(); inst_of <- c()
  for (ct in names(d)) {
    Y <- d[[ct]]$Y; lab <- rownames(Y); inst <- d[[ct]]$inst
    new <- lab
    for (i in unique(inst)) { idx <- which(inst == i); new[idx] <- lab[sample(idx)] }
    nc <- d[[ct]]$ncell
    rownames(Y) <- new; names(nc) <- new
    pb[[ct]] <- Y; ncell[[ct]] <- nc
    inst_of[new] <- as.character(inst)
  }
  u <- sort(unique(unlist(lapply(pb, rownames))))
  zs <- setNames(rnorm(length(u)), u)
  z_truth <- list()
  if (!is.null(pair)) {
    if (scenario == "S7") {
      g <- sample(Reduce(intersect, lapply(pair, function(ct) colnames(pb[[ct]]))), n_prog_genes)
      w_unit <- rnorm(n_prog_genes)
    }
    for (ct in pair) {
      Y <- pb[[ct]]
      z <- if (scenario == "S7") setNames(rnorm(length(u)), u) else zs
      if (scenario == "S7") { gc <- g; w <- w_unit * amp * d[[ct]]$msd }
      else { gc <- sample(colnames(Y), n_prog_genes); w <- rnorm(n_prog_genes, sd = amp * d[[ct]]$msd) }
      Y[, gc] <- Y[, gc] + outer(z[rownames(Y)], w)
      pb[[ct]] <- Y; z_truth[[ct]] <- z[rownames(Y)]
    }
  }
  meta <- data.frame(observation_id = u, institute = factor(inst_of[u]), stringsAsFactors = FALSE)
  list(x = as_cell_program_data(pb, sample_metadata = meta), ncell = ncell,
       truth = list(scenario = scenario, pair = pair, z = z_truth),
       strata = setNames(as.character(inst_of[u]), u))
}
