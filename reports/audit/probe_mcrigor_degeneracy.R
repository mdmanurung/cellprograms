# Probe: which membership degeneracy trips mcRigor's "No cell overlap" error?
# Cases mimic the integration fixture (60 genes x 50 cells per group) and a
# known-good control. Tracebacks captured on error.
set.seed(1)
mk_counts <- function(ng, nc) {
  m <- Matrix::Matrix(rpois(ng * nc, 2), nrow = ng, sparse = TRUE,
                      dimnames = list(paste0("gene", 1:ng), paste0("cell", 1:nc)))
}
mk_memb <- function(nc, sizes) {
  stopifnot(sum(sizes) == nc)
  rep(seq_along(sizes), sizes) |> as.character() |> setNames(paste0("cell", 1:nc))
}
probe <- function(label, counts, memb) {
  cat("\n==== ", label, " ====\n", sep = "")
  res <- tryCatch(
    cellprograms:::.mcrigor_detect_group(counts, memb, mode = "average",
                                         mcrigor_args = list(), group_id = "probe"),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); print(sys.calls()) ; NULL },
    warning = function(w) { cat("WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }
  )
  if (is.null(res)) cat("-> failed\n") else cat("-> OK, rows:", nrow(res), "\n")
  invisible(res)
}
pkgload::load_all("/workspace/repo", quiet = TRUE)

c50 <- mk_counts(60, 50)
probe("A: 2 balanced metacells (25/25)", c50, mk_memb(50, c(25, 25)))
probe("B: 7 metacells, 1 dominant (44/1/1/1/1/1/1)", c50, mk_memb(50, c(44, 1, 1, 1, 1, 1, 1)))
probe("C: single metacell (50)", c50, mk_memb(50, c(50)))
c120 <- mk_counts(200, 120)
probe("D: control 6x20 balanced", c120, mk_memb(120, rep(20, 6)))
