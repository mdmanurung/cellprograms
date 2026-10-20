# Probe: does mcRigor_OPTIMIZE hit the single-metacell edge case, and how does
# OptimizeGammaSEACells' fallback degrade?
# Fixture-like regime: 50-cell group, gammas 25/50/100 -> n_seacells = 2/1/1,
# i.e. one usable column and two single-valued columns.
pkgload::load_all("/workspace/repo", quiet = TRUE)
set.seed(1)
nc <- 50
counts <- Matrix::Matrix(rpois(60 * nc, 2), 60, nc, sparse = TRUE,
                         dimnames = list(paste0("gene", 1:60), paste0("cell", 1:nc)))
obj <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA",
                                  min.features = 0, min.cells = 0)
mk_col <- function(sizes) rep(seq_along(sizes), sizes) |> as.character()
# fixture regime: gamma 25 -> 2 metacells; gammas 50/100 -> 1 metacell
mem_bad <- data.frame(`25` = mk_col(c(25, 25)), `50` = mk_col(c(50)),
                      `100` = mk_col(c(50)), check.names = FALSE,
                      stringsAsFactors = FALSE)
rownames(mem_bad) <- paste0("cell", 1:nc)
# control: all columns >= 2 metacells
mem_ok <- data.frame(`25` = mk_col(c(13, 12, 13, 12)), `50` = mk_col(c(25, 25)),
                     `100` = mk_col(c(17, 17, 16)), check.names = FALSE,
                     stringsAsFactors = FALSE)
rownames(mem_ok) <- paste0("cell", 1:nc)

run <- function(label, mem) {
  cat("\n====", label, "====\n")
  t0 <- Sys.time()
  res <- tryCatch(
    mcRigor::mcRigor_OPTIMIZE(obj, cell_membership = mem, Nrep = 5,
                              gene_filter = 0.1, feature_use = 2000,
                              test_cutoff = 0.01, cor_method = "pearson",
                              aggregate_method = "mean",
                              fields = character(0),
                              draw = FALSE, step_save = FALSE,
                              output_file = tempfile(fileext = ".rds")),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("WARN:", conditionMessage(w), "\n"); NULL }
  )
  cat(sprintf("-> %s (%.1fs)\n", if (is.null(res)) "FAILED" else "OK",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(res)
}
r_bad <- run("fixture-like: gammas 25/50/100 = 2/1/1 metacells", mem_bad)
r_ok <- run("control: all gammas >= 2 metacells", mem_ok)
if (!is.null(r_ok)) {
  pick <- cellprograms:::.mcrigor_optimize_pick(r_ok, c(25, 50, 100))
  cat("control pick: gamma =", pick$gamma, " scores =", pick$scores, "\n")
}
