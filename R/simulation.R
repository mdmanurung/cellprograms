# cellprograms: simulation of pseudobulk data with known program structure
#
# Implements the ground-truth scenarios from the benchmarking plan
# (Sections 10-11). Used by the feasibility spike, the TusoAI evaluator,
# and package validation.

#' Simulate cell-type pseudobulk matrices with planted programs
#'
#' @param scenario One of "private" (Scenario 1), "mixed" (Scenario 4,
#'   principal), "same_genes_independent" (Scenario 7 negative control),
#'   or "global" (Scenario 2).
#' @param n_samples Number of biological observations.
#' @param n_genes Total gene pool size.
#' @param cell_types Character vector of cell-type names.
#' @param program_size Number of active genes per program.
#' @param signal_sd Program activity standard deviation.
#' @param noise_sd Residual noise standard deviation.
#' @param seed Random seed.
#'
#' @return A list with `pseudobulk` (named list of observations x genes
#'   matrices), `truth` (per-cell-type list of true score vectors, loading
#'   indices, and program labels), and `scenario`.
#'
#' @export
simulate_pseudobulk <- function(scenario = "mixed", n_samples = 150,
                                n_genes = 1500, program_size = 50,
                                signal_sd = 1.5, noise_sd = 1.0,
                                cell_types = c("B", "CD4", "CD8", "NK", "Mono"),
                                seed = 1) {
  set.seed(seed)
  stopifnot(scenario %in% c("private", "mixed", "same_genes_independent", "global"))

  gene_pool <- paste0("G", seq_len(n_genes))
  truth <- list()
  pseudobulk <- list()

  # helper: add a program to a cell type's expression
  add_program <- function(Y, label, gene_idx, z) {
    w <- rep(0, ncol(Y))
    w[gene_idx] <- rnorm(length(gene_idx), 0, 1)
    Y <- Y + z %*% matrix(w, nrow = 1)
    list(Y = Y, label = label, gene_idx = gene_idx, z = z, w = w)
  }

  for (ct in cell_types) {
    Y <- matrix(0, nrow = n_samples, ncol = n_genes,
                dimnames = list(paste0("obs", seq_len(n_samples)), gene_pool))
    truth[[ct]] <- list()
    prog_id <- 0
    add <- function(label, gene_idx, z) {
      prog_id <<- prog_id + 1
      res <- add_program(Y, label, gene_idx, z)
      Y <<- res$Y
      truth[[ct]][[prog_id]] <<- list(
        label = label, gene_idx = gene_idx, z = res$z, w = res$w
      )
    }

    if (scenario == "private") {
      g <- sample(n_genes, program_size)
      add(paste0(ct, "_P1"), g, rnorm(n_samples, 0, signal_sd))
      g <- sample(n_genes, program_size)
      add(paste0(ct, "_P2"), g, rnorm(n_samples, 0, signal_sd))
    } else if (scenario == "mixed") {
      # B-private, Mono-private, CD8+NK shared, global IFN
      if (ct == "B") {
        add("B_private", sample(n_genes, program_size), rnorm(n_samples, 0, signal_sd))
      }
      if (ct == "Mono") {
        add("Mono_private", sample(n_genes, program_size), rnorm(n_samples, 0, signal_sd))
      }
      if (ct %in% c("CD8", "NK")) {
        # shared program: same activity vector, cell-type-specific genes
        if (ct == "CD8") {
          shared_z <- rnorm(n_samples, 0, signal_sd)
          shared_g_cd8 <- sample(n_genes, program_size)
        }
        add("CD8NK_shared", shared_g_cd8, shared_z)
      }
      if (ct %in% c("B", "CD4", "NK", "Mono")) {
        # one shared activity across its member cell types, cell-type-specific genes
        if (ct == "B") ifn_z <- rnorm(n_samples, 0, signal_sd)
        add("global_IFN", sample(n_genes, program_size), ifn_z)
      }
    } else if (scenario == "same_genes_independent") {
      # B and Mono share the SAME gene set (IFN-like) but have INDEPENDENT
      # sample activities; other cell types get unrelated programs.
      if (ct %in% c("B", "Mono")) {
        if (ct == "B") ifn_genes <- sample(n_genes, program_size)
        add("IFN_like", ifn_genes, rnorm(n_samples, 0, signal_sd))
      } else {
        add(paste0(ct, "_P1"), sample(n_genes, program_size), rnorm(n_samples, 0, signal_sd))
      }
    } else if (scenario == "global") {
      if (ct == cell_types[1]) global_z <- rnorm(n_samples, 0, signal_sd)
      add("global_state", sample(n_genes, program_size), global_z)
    }

    # residual noise
    Y <- Y + matrix(rnorm(n_samples * n_genes, 0, noise_sd), n_samples, n_genes)
    pseudobulk[[ct]] <- Y
  }

  list(pseudobulk = pseudobulk, truth = truth, scenario = scenario,
       n_samples = n_samples, n_genes = n_genes, seed = seed)
}
