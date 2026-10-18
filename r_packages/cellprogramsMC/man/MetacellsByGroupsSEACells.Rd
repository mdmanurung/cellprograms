\name{MetacellsByGroupsSEACells}
\alias{MetacellsByGroupsSEACells}
\title{SEACells-backed metacells per group with mcRigor QC}
\usage{
MetacellsByGroupsSEACells(
  seurat_obj, group.by = "seurat_clusters", ident.group = NULL,
  reduction = "pca", dims = NULL, assay = "RNA", slot = "counts",
  layer = "counts", mode = c("sum", "average"), cells_per_metacell = 75,
  target_metacells = NULL, min_cells = 100, n_top_genes = 1500, n_pc = 50,
  n_waypoint_eigs = 10, waypt_proportion = 0.9, n_neighbors = 15,
  convergence_epsilon = 1e-5, max_iter = 20, use_sparse = TRUE,
  use_soft_assignments = TRUE, run_mcrigor = TRUE, mcrigor_args = list(),
  seed = 20261016, name = "SEAC", wgcna_name = "SEACells", verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object.}
\item{group.by}{Metadata column(s) defining the groups (e.g. donor, cell
type, condition). One SEACells model is fit per group.}
\item{ident.group}{Optional prefix for metacell ids; defaults to the
\code{group.by} columns joined by underscores.}
\item{reduction}{Dimensional reduction to use for kernel construction;
\code{"pca"} by default. Set to \code{NULL} to build the kernel on the
expression matrix directly.}
\item{dims}{Optional PCA dimensions to keep (default: all available).}
\item{assay, slot, layer}{Assay and counts layer (Seurat v5 layers or v4
slots are both supported).}
\item{mode}{Aggregation mode for metacell profiles: \code{"sum"} (default;
sums counts within a metacell) or \code{"average"}. Note the aggregation-
scale caveat: sum-aggregated profiles live on a different nonlinear scale
after \code{log1p} than per-cell profiles (M6/E6 G4).}
\item{cells_per_metacell}{Target granularity gamma: approximate number of
cells per metacell.}
\item{target_metacells}{Optional explicit number of metacells per group
(overrides \code{cells_per_metacell}).}
\item{min_cells}{Groups with fewer cells are skipped (recorded in
\code{skipped_groups}).}
\item{n_top_genes}{Number of highly variable genes used by SEACells.}
\item{n_pc}{Number of PCA dimensions used by SEACells.}
\item{n_waypoint_eigs, waypt_proportion, n_neighbors}{SEACells kernel
parameters.}
\item{convergence_epsilon, max_iter}{SEACells optimization stopping controls.}
\item{use_sparse}{Kept for API compatibility; SEACells 0.3.3 removed the
constructor kwarg.}
\item{use_soft_assignments}{Use the soft-assignment matrix for aggregation
(hard labels are reconstructed for QC).}
\item{run_mcrigor}{Run mcRigor QC per group (mcDiv, dubious verdicts,
ZeroRate/DubRate/Score).}
\item{mcrigor_args}{Optional named list overriding mcRigor arguments
(\code{Nrep}, \code{gene_filter}, \code{feature_use}, \code{test_cutoff},
\code{cor_method}, ...).}
\item{seed}{Seed forwarded to the Python side (\code{np.random.RandomState}).}
\item{name}{Prefix for metacell ids.}
\item{wgcna_name}{misc slot name under which results are stored.}
\item{verbose}{Print progress messages.}
}
\value{
The Seurat object with results stored under
\code{misc[[wgcna_name]]}: \code{metacell_object} (genes x metacells
aggregated profiles), \code{soft_assignments}, \code{qc} (per-metacell
mcRigor table), \code{fit_stats}, \code{zerorate_overall},
\code{dubrate_overall}, \code{score_overall}, \code{params},
\code{skipped_groups}, \code{call}.
}
\description{
Constructs transcriptomically homogeneous metacells from a Seurat object
using SEACells (Persad et al. 2023) as the aggregation backend, stratified
by user-supplied groupings following the hdWGCNA \code{MetacellsByGroups}
interface (Morabito et al. 2023). Metacell quality is governed by mcRigor
(Liu & Li 2025): dubious-metacell detection via mcDiv double-permutation
nulls with size-conditional thresholds, plus DubRate/ZeroRate scoring.
}
\details{
SEACells is GPLv2 and mcRigor is MIT; both run as external dependencies
inside an isolated basilisk environment and are never source-copied into
this package.

Transfer caveats (M6/E6, registered gates): metacells trade rare-state
fidelity (G1) and soft-projection exactness (G4) for denoising and scale.
Downstream program discovery on the metacell route should use a fixed K
with pruning and report ZeroRate and effective K per fit (M6 revision R3);
never present metacell-route rare-program results without the transfer
caveat.
}
\examples{
\dontrun{
obj <- MetacellsByGroupsSEACells(obj, group.by = c("donor", "cell_type"),
                                 cells_per_metacell = 75)
}
}
