## simulate.R — ground-truth simulation suite (benchmarking plan Scenarios 1-7)
##
## Each scenario generates sample x gene x cell-type pseudobulk-like matrices
## with known program structure: private, partially shared, global,
## correlated-but-distinct, and negative-control configurations.

#' Simulate a cell-type x sample x gene dataset with known programs
#'
#' @param scenario Integer 1-7 (see package benchmarking plan).
#' @param N Number of observations (samples).
#' @param G Genes per cell type.
#' @param cell_types Cell-type labels.
#' @param signal_sd SD of active gene loadings.
#' @param sparsity Fraction of genes active in each program.
#' @param noise_sd Noise SD (per entry).
#' @param cor_strength Activity correlation for Scenario 5.
#' @param missing_frac Fraction of sample x cell-type combinations to drop.
#' @param seed Random seed.
#' @return A list with `x` (cell_program_data) and `truth` (true activities,
#'   loadings, and program structure table).
#' @export
simulate_programs <- function(scenario = 4L, N = 100L, G = 500L,
                              cell_types = c("B", "CD4", "CD8", "NK", "Mono"),
                              signal_sd = 3, sparsity = 0.05, noise_sd = 1,
                              cor_strength = 0.7, missing_frac = 0,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- length(cell_types)
  obs_ids <- paste0("obs", seq_len(N))
  genes <- paste0("g", seq_len(G))

  ## helper: sparse loading vector
  new_loadings <- function(active_genes = NULL) {
    w <- rep(0, G)
    idx <- active_genes %||% sample.int(G, max(2L, round(G * sparsity)))
    w[idx] <- stats::rnorm(length(idx), 0, signal_sd)
    stats::setNames(w, genes)
  }

  activities <- list()  # name -> list(z, cell_types, type)
  loadings <- stats::setNames(vector("list", C), cell_types)
  for (ct in cell_types) loadings[[ct]] <- list()

  add_program <- function(name, cts, type, z = NULL, shared_genes = NULL) {
    if (is.null(z)) z <- stats::rnorm(N)
    activities[[name]] <<- list(z = z, cell_types = cts, type = type)
    for (ct in cts) {
      loadings[[ct]][[name]] <<- new_loadings(shared_genes)
    }
  }

  switch(
    as.character(scenario),
    "1" = {
      ## purely private programs, independent activities
      for (ct in cell_types) {
        add_program(paste0(ct, "_P1"), ct, "private")
        add_program(paste0(ct, "_P2"), ct, "private")
      }
    },
    "2" = {
      ## one global program (same activity, CT-specific genes) + one private each
      zg <- stats::rnorm(N)
      add_program("GLOBAL_1", cell_types, "global", z = zg)
      for (ct in cell_types) add_program(paste0(ct, "_P1"), ct, "private")
    },
    "3" = {
      ## lineage-restricted sharing: CD8 + NK
      zs <- stats::rnorm(N)
      add_program("CD8NK_shared", c("CD8", "NK"), "partial", z = zs)
      for (ct in cell_types) add_program(paste0(ct, "_P1"), ct, "private")
    },
    "4" = {
      ## mixed: B-private, Mono-private, CD8/NK-shared, global IFN-like
      add_program("B_priv", "B", "private")
      add_program("Mono_priv", "Mono", "private")
      add_program("CD8NK_shared", c("CD8", "NK"), "partial")
      add_program("IFN_global", cell_types, "global")
    },
    "5" = {
      ## correlated but distinct programs in B and Mono
      z1 <- stats::rnorm(N)
      z2 <- cor_strength * z1 + sqrt(1 - cor_strength^2) * stats::rnorm(N)
      add_program("B_distinct", "B", "private", z = z1)
      add_program("Mono_distinct", "Mono", "private", z = z2)
      for (ct in setdiff(cell_types, c("B", "Mono"))) {
        add_program(paste0(ct, "_P1"), ct, "private")
      }
    },
    "6" = {
      ## same systemic state, entirely different genes per cell type
      z <- stats::rnorm(N)
      add_program("SYS", cell_types, "global", z = z)
      for (ct in cell_types) add_program(paste0(ct, "_P1"), ct, "private")
    },
    "7" = {
      ## negative control: same genes in B and Mono, independent activities
      shared_idx <- sample.int(G, max(2L, round(G * sparsity)))
      add_program("B_ifn_like", "B", "private", shared_genes = shared_idx)
      add_program("Mono_ifn_like", "Mono", "private", shared_genes = shared_idx)
      for (ct in setdiff(cell_types, c("B", "Mono"))) {
        add_program(paste0(ct, "_P1"), ct, "private")
      }
    },
    .stopf("Unknown scenario %d.", scenario)
  )

  ## build expression matrices
  mats <- stats::setNames(vector("list", C), cell_types)
  for (ct in cell_types) {
    Y <- matrix(0, N, G, dimnames = list(obs_ids, genes))
    for (pn in names(loadings[[ct]])) {
      z <- activities[[pn]]$z
      w <- loadings[[ct]][[pn]]
      Y <- Y + z %*% t(w)
    }
    Y <- Y + matrix(stats::rnorm(N * G, 0, noise_sd), N, G)
    mats[[ct]] <- Y
  }

  ## missingness: drop random sample x cell-type combinations
  if (missing_frac > 0) {
    for (ct in cell_types) {
      drop <- sample.int(N, floor(N * missing_frac))
      if (length(drop) > 0L) mats[[ct]] <- mats[[ct]][-drop, , drop = FALSE]
    }
  }

  truth <- list(
    activities = activities,
    loadings = loadings,
    structure = do.call(rbind, lapply(names(activities), function(pn) {
      data.frame(program = pn, type = activities[[pn]]$type,
                 cell_types = paste(activities[[pn]]$cell_types, collapse = "+"),
                 stringsAsFactors = FALSE)
    }))
  )

  list(x = as_cell_program_data(mats), truth = truth)
}
