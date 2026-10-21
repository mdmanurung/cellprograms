# Stage 1: synthetic fixtures with planted structure + injected doublets.
# Usage: Rscript benchmark_fixture.R <out_dir> <scale: 5k|15k> <seed>
# Mirrors the mcRigor_reproduce synthetic design (planted clusters, injected
# doublets with ground truth, donor structure) at package-relevant scale.

suppressPackageStartupMessages({
  library(Matrix)
})

build_fixture <- function(n_per_type, n_types = 5L, n_donors = 10L, n_genes = 3000L,
                          program_size = 50L, doublet_rate = 0.05, seed = 1) {
  set.seed(seed)
  cell_types <- paste0("T", seq_len(n_types))
  donors <- paste0("D", seq_len(n_donors))
  genes <- paste0("gene_", seq_len(n_genes))

  # Type base profiles: sparse random log-mean profiles (most genes off-ish per type)
  base_mu <- matrix(0, nrow = n_genes, ncol = n_types, dimnames = list(genes, cell_types))
  base_mu[,] <- exp(rnorm(n_genes * n_types, mean = 0.2, sd = 1.0))
  # make profiles sparse-ish and type-discriminating
  type_specific <- matrix(rbinom(n_genes * n_types, 1, 0.15), n_genes, n_types) == 1
  base_mu[type_specific] <- base_mu[type_specific] * exp(rnorm(sum(type_specific), 1.2, 0.5))
  base_mu[!type_specific] <- base_mu[!type_specific] * 0.15

  # Donor effects (global shift per donor, small)
  donor_shift <- rnorm(n_donors, 0, 0.15)
  names(donor_shift) <- donors

  # Planted programs: per cell type, `program_size` genes driven by a shared
  # donor-level activity (the recovery target downstream)
  program_genes <- lapply(cell_types, function(ct) sample(genes, program_size))
  names(program_genes) <- cell_types
  donor_activity <- lapply(cell_types, function(ct) setNames(rnorm(n_donors, 0, 1), donors))
  names(donor_activity) <- cell_types

  sim_one <- function(ct) {
    mu <- base_mu[, ct]
    prog_g <- program_genes[[ct]]
    act <- donor_activity[[ct]]
    cell_donors <- sample(donors, n_per_type, replace = TRUE)
    cells <- lapply(seq_len(n_per_type), function(i) {
      d <- cell_donors[i]
      lib <- exp(rnorm(1, log(5000), 0.3))
      lam <- mu * lib * exp(donor_shift[d]) / sum(mu)
      # planted program: donor activity scales the program genes
      lam[prog_g] <- lam[prog_g] * exp(0.8 * act[d])
      rnbinom(n_genes, size = 2, mu = lam)
    })
    list(counts = do.call(cbind, cells), donors = cell_donors)
  }

  sims <- lapply(cell_types, sim_one)
  counts <- do.call(cbind, lapply(sims, `[[`, "counts"))
  donors_vec <- unlist(lapply(sims, `[[`, "donors"), use.names = FALSE)
  colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  rownames(counts) <- genes
  meta <- data.frame(
    cell_id = colnames(counts),
    cell_type = rep(cell_types, each = n_per_type),
    donor = donors_vec,
    stringsAsFactors = FALSE
  )

  # Inject doublets: replace 5% of cells with sums of two cells from DIFFERENT types
  n_doublets <- round(doublet_rate * ncol(counts))
  set.seed(seed + 2L)
  dbl_idx <- sample(ncol(counts), n_doublets)
  origins <- character(length(dbl_idx))
  for (k in seq_along(dbl_idx)) {
    j <- dbl_idx[k]
    types <- sample(cell_types, 2)
    c1 <- sample(colnames(counts)[meta$cell_type == types[1]], 1)
    c2 <- sample(colnames(counts)[meta$cell_type == types[2]], 1)
    origins[k] <- paste(sort(c(types[1], types[2])), collapse = "/")
    counts[, j] <- counts[, c1] + counts[, c2]
  }
  meta$dubious <- FALSE
  meta$doublet_origin <- NA_character_
  meta$dubious[dbl_idx] <- TRUE
  meta$doublet_origin[dbl_idx] <- origins

  list(counts = Matrix(counts, sparse = TRUE), meta = meta,
       program_genes = program_genes, donor_activity = donor_activity,
       genes = genes, cell_types = cell_types, donors = donors)
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; scale <- args[2]; seed <- as.integer(args[3])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_per_type <- switch(scale, "5k" = 1000L, "15k" = 3000L, stop("scale must be 5k or 15k"))
fx <- build_fixture(n_per_type = n_per_type, seed = seed)

saveRDS(fx$counts, file.path(out_dir, paste0("counts_", scale, "_s", seed, ".rds")))
saveRDS(fx$meta, file.path(out_dir, paste0("meta_", scale, "_s", seed, ".rds")))
saveRDS(fx$program_genes, file.path(out_dir, paste0("programs_", scale, "_s", seed, ".rds")))
saveRDS(fx$donor_activity, file.path(out_dir, paste0("activity_", scale, "_s", seed, ".rds")))
cat(sprintf("FIXTURE OK %s seed=%d: %d cells x %d genes, %d doublets\n",
            scale, seed, ncol(fx$counts), nrow(fx$counts), sum(fx$meta$dubious)))
