# cellprogramsMC package
#
# Basilisk environment for the SEACells backend.
#
# LICENSE NOTE (binding): SEACells (dpeerlab/SEACells) is GPLv2 per its
# LICENSE.txt (the PyPI trove classifier claiming MIT is contradicted by the
# license file and is treated as erroneous). It is therefore used strictly as
# an external dependency inside this isolated basilisk environment, never
# source-copied into this package. mcRigor (JSB-UCLA/mcRigor) is MIT and is a
# direct R dependency.

#' Basilisk environment for SEACells
#'
#' The Python environment used by [MetacellsByGroupsSEACells()] to run the
#' SEACells aggregation backend. The environment is created automatically on
#' first use (requires network access at that point).
#'
#' SEACells 0.3.3 (PyPI) requires Python >= 3.11; its core dependencies are
#' CPU-only (numpy, pandas, scipy, scanpy, anndata, palantir, numba,
#' scikit-learn, pyranges, matplotlib, seaborn, tqdm) and are resolved by pip.
#'
#' @format A `basilisk::BasiliskEnvironment` object.
"seacells_env"

seacells_env <- basilisk::BasiliskEnvironment(
    envname = "seacells_env",
    pkgname = "cellprogramsMC",
    packages = c("python=3.11"),
    # ipywidgets: SEACells imports tqdm.notebook (progress bars), which raises
    # ImportError("IProgress not found") in headless envs without it.
    pip = c("SEACells==0.3.3", "ipywidgets==8.1.7")
)

#' cellprogramsMC: SEACells-backed metacells with mcRigor QC
#'
#' Constructs metacells from a Seurat object with SEACells as the aggregation
#' backend, stratified by groups (donor, cell type, condition) following the
#' hdWGCNA `MetacellsByGroups` interface, and governs metacell quality with
#' mcRigor (mcDiv dubious-metacell detection, DubRate/ZeroRate scoring,
#' optional granularity optimization). Results are stored on the Seurat object
#' and can feed cellprograms / cebmf_torch program discovery.
#'
#' @docType package
#' @name cellprogramsMC-package
NULL
