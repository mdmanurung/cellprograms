\name{OptimizeGammaSEACells}
\alias{OptimizeGammaSEACells}
\title{Optimize SEACells granularity with mcRigor}
\usage{
OptimizeGammaSEACells(
  seurat_obj, group.by = "seurat_clusters", gammas = c(25, 50, 100),
  reduction = "pca", dims = NULL, assay = "RNA", slot = "counts",
  layer = "counts", n_top_genes = 1500, n_pc = 50, n_waypoint_eigs = 10,
  waypt_proportion = 0.9, n_neighbors = 15, convergence_epsilon = 1e-5,
  max_iter = 20, use_sparse = TRUE, min_cells = 100, mcrigor_args = list(),
  seed = 20261016, verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object.}
\item{group.by}{Metadata column(s) defining groups (as in
\code{\link{MetacellsByGroupsSEACells}}).}
\item{gammas}{Candidate granularities as cells-per-metacell values. Default
\code{c(25, 50, 100)} (the mcRigor semi-synthetic optimum was 50 for
SEACells).}
\item{reduction, dims, assay, slot, layer, n_top_genes, n_pc,
n_waypoint_eigs, waypt_proportion, n_neighbors, convergence_epsilon,
max_iter, use_sparse, min_cells}{Passed to the
\code{\link{MetacellsByGroupsSEACells}} internals.}
\item{mcrigor_args}{Optional named list overriding mcRigor_OPTIMIZE
arguments (\code{D_bw}, \code{optim_method}, \code{weight}, \code{Nrep}, ...).}
\item{seed}{Seed forwarded to the Python side.}
\item{verbose}{Print progress messages.}
}
\value{
A list with \code{optimal_gamma} (per group), \code{scores} (per group x
gamma Score table), and \code{membership} (the per-cell membership
data.frame for the chosen partition, usable with
\code{mcRigor::mcRigorTS_Step1}/\code{mcRigorTS_Step2} for two-step
re-partitioning of dubious metacells).
}
\description{
Runs SEACells at several granularities (cells per metacell) per group and
calls \code{mcRigor::mcRigor_OPTIMIZE} to select the granularity maximizing
the mcRigor Score (1 - w*DubRate - (1-w)*ZeroRate, w = 0.5 default).
Expensive: one SEACells fit per group per candidate granularity.
}
\examples{
\dontrun{
opt <- OptimizeGammaSEACells(obj, group.by = c("donor", "cell_type"),
                             gammas = c(25, 50, 100))
}
}
