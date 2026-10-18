# Basilisk environment for the SEACells metacell backend.
#
# Folded in from the retired cellprogramsMC sub-package (ledger 2026-10-18).
# Follows the same lazy pattern as .cebmf_basilisk_env() in ebmf_backend.R:
# basilisk stays in Suggests and the environment object is created on first
# use, so installing cellprograms never requires basilisk or Python.
#
# LICENSE NOTE (binding): SEACells (dpeerlab/SEACells) is GPLv2 per its
# LICENSE.txt (the PyPI trove classifier claiming MIT is contradicted by the
# license file and is treated as erroneous). It is therefore used strictly as
# an external dependency inside this isolated basilisk environment, never
# source-copied into this package. mcRigor (JSB-UCLA/mcRigor) is MIT and is a
# suggested R dependency.

#' Basilisk environment for the SEACells backend
#'
#' The Python environment used by [MetacellsByGroupsSEACells()] to run the
#' SEACells aggregation backend. The environment is created automatically on
#' first use (requires network access at that point).
#'
#' SEACells 0.3.3 (PyPI) requires Python >= 3.11; its core dependencies are
#' CPU-only (numpy, pandas, scipy, scanpy, anndata, palantir, numba,
#' scikit-learn, pyranges, matplotlib, seaborn, tqdm) and are resolved by pip.
#'
#' @return A `basilisk::BasiliskEnvironment` object.
#' @noRd
.seacells_env <- function() {
  if (!requireNamespace("basilisk", quietly = TRUE)) {
    cli::cli_abort(
      "{.pkg basilisk} is required for SEACells metacell construction. Install it with {.code BiocManager::install('basilisk')}."
    )
  }
  basilisk::BasiliskEnvironment(
    envname = "seacells_env",
    pkgname = "cellprograms",
    packages = c("python=3.11"),
    # ipywidgets: SEACells imports tqdm.notebook (progress bars), which raises
    # ImportError("IProgress not found") in headless envs without it.
    pip = c("SEACells==0.3.3", "ipywidgets==8.1.7")
  )
}
