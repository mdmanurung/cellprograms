"""PATpy adapter for the cellprograms EBMF sample representation.

Exposes per-cell-type EBMF program scores (flashier backend, canonicalized)
as a patient/sample representation compatible with patpy's
``SampleRepresentationMethod`` interface, inheriting patpy's published metric
suite (kNN, distances-significance, silhouette, persistence, PERMANOVA,
linear probe).

The representation is the column-concatenation of per-cell-type canonical
score matrices (samples x K_ct). Each cell-type block is z-scaled per column
(``block_scale="z"``) so cell types with many factors do not dominate the
euclidean distance. Samples missing a cell type receive an all-zero block
(blocks are mean-zero by construction). Cell types where EBMF retained zero
factors contribute no columns.

Requires an R installation with the ``cellprograms`` package (and its
flashier/ebnm dependencies) available. Set ``R_HOME`` and ``R_LIBS_USER``
before importing if rpy2 is not yet loaded.

Example
-------
>>> from patpy_cellprograms import EBMF
>>> method = EBMF(sample_key="sample", cell_group_key="cell_type")
>>> method.prepare_anndata(adata)
>>> distances = method.calculate_distance_matrix()
>>> from patpy.tl import evaluate_representation
>>> evaluate_representation(distances, adata.obs.groupby("sample")["disease"].first(), method="knn")
"""

from __future__ import annotations

import os
import warnings
from uuid import uuid4

import numpy as np
import pandas as pd
import scipy.spatial

from patpy.tl.sample_representation import SampleRepresentationMethod, valid_distance_metric

DEFAULT_R_LIBS = os.environ.get("R_LIBS_USER", "/workspace/.Rlib")


def _ensure_rpy2():
    """Import rpy2, defaulting R_HOME/R_LIBS_USER to the conda R + cellprograms lib."""
    if "rpy2" in __import__("sys").modules:
        return __import__("rpy2.robjects", fromlist=["robjects"])
    os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
    os.environ.setdefault("R_LIBS_USER", DEFAULT_R_LIBS)
    import rpy2.robjects as ro

    return ro


class EBMF(SampleRepresentationMethod):
    """Sample representation via cell-type-specific EBMF (cellprograms R package).

    Parameters
    ----------
    sample_key : str
        Column in ``adata.obs`` containing sample (donor) IDs.
    cell_group_key : str
        Column in ``adata.obs`` containing cell-type labels.
    layer : str or None
        Feature source passed to the base class (``None`` → ``adata.X``).
    seed : int
        Random seed forwarded to ``fit_celltype_programs``.
    loading_prior : str
        Gene-side EBNM prior family: ``"point_laplace"`` (default),
        ``"point_normal"``, or ``"unimodal"``.
    max_factors : int
        Maximum candidate factors per cell type (greedy_Kmax).
    var_type : int
        flashier variance type (locked default 1).
    backfit, nullcheck : bool
        flashier refinement options (locked defaults True).
    features : {"all", "variable"}
        Feature selection mode (locked default ``"all"``).
    n_variable_genes : int
        Number of most-variable genes when ``features="variable"``.
    block_scale : {"z", "none"}
        Per-cell-type block scaling of score columns before concatenation.
    dist : str
        Distance metric for the sample-sample matrix (default ``"euclidean"``).
    r_lib_paths : str or None
        Extra R library path containing ``cellprograms`` (default
        ``R_LIBS_USER`` or ``/workspace/.Rlib``).
    """

    def __init__(
        self,
        sample_key: str,
        cell_group_key: str,
        layer: str | None = None,
        seed: int = 67,
        loading_prior: str = "point_laplace",
        max_factors: int = 30,
        var_type: int = 1,
        backfit: bool = True,
        nullcheck: bool = True,
        features: str = "all",
        n_variable_genes: int = 2000,
        block_scale: str = "z",
        dist: str = "euclidean",
        r_lib_paths: str | None = None,
    ):
        super().__init__(sample_key=sample_key, cell_group_key=cell_group_key, layer=layer, seed=seed)
        if block_scale not in ("z", "none"):
            raise ValueError(f"block_scale must be 'z' or 'none', got {block_scale!r}")
        self.loading_prior = loading_prior
        self.max_factors = max_factors
        self.var_type = var_type
        self.backfit = backfit
        self.nullcheck = nullcheck
        self.features = features
        self.n_variable_genes = n_variable_genes
        self.block_scale = block_scale
        self.dist = dist
        self.r_lib_paths = r_lib_paths
        self.sample_representation: pd.DataFrame | None = None
        self._r_fit_id: str | None = None

    # ------------------------------------------------------------------ fit
    def prepare_anndata(self, adata) -> None:
        """Pseudobulk per (sample, cell type) and fit EBMF in R."""
        super().prepare_anndata(adata=adata)
        if self.cell_group_key is None or self.cell_group_key not in adata.obs.columns:
            raise ValueError(f"cell_group_key='{self.cell_group_key}' is required for EBMF.")

        X = self._get_data()
        sparse = hasattr(X, "tocsr")
        obs = pd.DataFrame(
            {
                "sample": adata.obs[self.sample_key].astype(str).values,
                "cell_type": adata.obs[self.cell_group_key].astype(str).values,
            },
            index=adata.obs_names.astype(str),
        )
        genes = np.asarray(adata.var_names.astype(str))

        # Pseudobulk: mean expression per (sample, cell type). Sparse-aware:
        # indicator-matrix product avoids densifying the full cell matrix.
        sample_labels = sorted(obs["sample"].unique())
        sample_pos = {s: i for i, s in enumerate(sample_labels)}
        pseudobulk: dict[str, pd.DataFrame] = {}
        for ct in sorted(obs["cell_type"].unique()):
            mask = (obs["cell_type"] == ct).values
            sub = X[mask]
            rows = np.array([sample_pos[s] for s in obs["sample"].values[mask]])
            n_s = len(sample_labels)
            if sparse:
                from scipy import sparse as _sp

                A = _sp.csr_matrix(
                    (np.ones(mask.sum()), (rows, np.arange(mask.sum()))),
                    shape=(n_s, mask.sum()),
                )
                agg = A @ sub  # sums per sample
                counts = np.bincount(rows, minlength=n_s)
                agg = agg.toarray() / np.maximum(counts, 1)[:, None]
            else:
                agg = np.zeros((n_s, X.shape[1]))
                np.add.at(agg, rows, np.asarray(sub))
                counts = np.bincount(rows, minlength=n_s)
                agg = agg / np.maximum(counts, 1)[:, None]
            keep = counts > 0
            mat = pd.DataFrame(agg[keep], index=np.array(sample_labels)[keep], columns=genes)
            pseudobulk[ct] = mat

        self._fit_in_r(pseudobulk)
        self._fitted = True

    def _fit_in_r(self, pseudobulk: dict[str, pd.DataFrame]) -> None:
        ro = _ensure_rpy2()
        from rpy2.robjects import numpy2ri, pandas2ri

        lib_paths = self.r_lib_paths or os.environ.get("R_LIBS_USER", DEFAULT_R_LIBS)
        fit_id = f"cpd_fit_{uuid4().hex[:12]}"
        self._r_fit_id = fit_id
        r_backfit = "TRUE" if self.backfit else "FALSE"
        r_nullcheck = "TRUE" if self.nullcheck else "FALSE"

        with (ro.default_converter + numpy2ri.converter + pandas2ri.converter).context():
            ro.r(f'.libPaths(c("{lib_paths}", .libPaths())); suppressMessages(library(cellprograms))')
            ro.globalenv[f"{fit_id}_cts"] = ro.StrVector(list(pseudobulk.keys()))
            # Index-based variable names (1-aligned with R's seq_along below):
            # cell-type labels can contain hyphens or other characters that
            # are invalid in R symbols.
            for i, (ct, mat) in enumerate(pseudobulk.items()):
                ro.globalenv[f"{fit_id}_Y_{i + 1}"] = mat.values
                ro.globalenv[f"{fit_id}_obs_{i + 1}"] = ro.StrVector(mat.index.tolist())
                ro.globalenv[f"{fit_id}_genes_{i + 1}"] = ro.StrVector(mat.columns.tolist())
            ro.r(
                f"""
                {fit_id} <- list()
                for (i in seq_along({fit_id}_cts)) {{
                    Y <- get(paste0("{fit_id}_Y_", i))
                    dimnames(Y) <- list(
                        get(paste0("{fit_id}_obs_", i)),
                        get(paste0("{fit_id}_genes_", i))
                    )
                    {fit_id}[[i]] <- Y
                }}
                names({fit_id}) <- {fit_id}_cts
                {fit_id}_data <- as_cell_program_data({fit_id})
                {fit_id}_fit <- canonicalize_programs(fit_celltype_programs(
                    {fit_id}_data,
                    loading_prior = "{self.loading_prior}",
                    max_factors = {int(self.max_factors)},
                    var_type = {int(self.var_type)},
                    backfit = {r_backfit},
                    nullcheck = {r_nullcheck},
                    features = "{self.features}",
                    n_variable_genes = {int(self.n_variable_genes)},
                    seed = {int(self.seed)}
                ))
                """
            )
            # Positional [[...]] access: `$`-chaining breaks on cell-type
            # names that are not valid R symbols (e.g. "B_non-switched_memory"
            # parses `$B_non-switched_memory` as subtraction).
            self._scores = {
                ct: pd.DataFrame(np.asarray(ro.r(f'{fit_id}_fit[["scores"]][[{i + 1}]]')))
                for i, ct in enumerate(pseudobulk)
            }
            # Recover labels: R matrix dimnames survive conversion.
            self._score_index = {
                ct: [str(x) for x in ro.r(f'rownames({fit_id}_fit[["scores"]][[{i + 1}]])')]
                for i, ct in enumerate(pseudobulk)
            }
            self._score_cols = {
                ct: [str(x) for x in ro.r(f'colnames({fit_id}_fit[["scores"]][[{i + 1}]])')]
                for i, ct in enumerate(pseudobulk)
            }
            # Clean globalenv except the fit (kept for loadings access).
            ro.r(f"rm(list = ls()[grepl('{fit_id}_(Y|obs|genes|cts)', ls())])")

    # ------------------------------------------------------------ transform
    def build_representation(self) -> pd.DataFrame:
        """Concatenate z-scaled per-cell-type score blocks on the union of samples."""
        blocks = []
        for ct, scores in self._scores.items():
            K = scores.shape[1]
            if K == 0:
                continue  # K=0 guard: cell type retained no factors
            mat = pd.DataFrame(scores.values, index=self._score_index[ct], columns=self._score_cols[ct])
            if self.block_scale == "z":
                sd = mat.std(axis=0, ddof=0)
                sd = sd.replace(0.0, 1.0)
                mat = (mat - mat.mean(axis=0)) / sd
            blocks.append(mat)

        if not blocks:
            raise RuntimeError("EBMF retained zero factors in every cell type; no representation available.")

        all_samples = sorted({idx for block in blocks for idx in block.index})
        rep = pd.DataFrame(0.0, index=all_samples, columns=pd.Index([], dtype=str))
        parts = [block.reindex(all_samples).fillna(0.0) for block in blocks]
        rep = pd.concat(parts, axis=1)
        rep.columns = rep.columns.astype(str)
        return rep

    def calculate_distance_matrix(self, force: bool = False) -> np.ndarray:
        """Sample-sample distance matrix on the concatenated EBMF representation."""
        distances = super().calculate_distance_matrix(force=force)
        if distances is not None:
            return distances

        rep = self.build_representation()
        self.sample_representation = rep
        metric = valid_distance_metric(self.dist)
        if metric == "euclidean":
            distances = scipy.spatial.distance.squareform(
                scipy.spatial.distance.pdist(rep.values, metric="euclidean")
            )
        else:
            distances = scipy.spatial.distance.squareform(
                scipy.spatial.distance.pdist(rep.values, metric=metric)
            )
        self._distances = distances
        self.samples = np.array(rep.index)
        return distances

    def get_loadings(self) -> dict[str, pd.DataFrame]:
        """Return per-cell-type canonical gene loadings (genes x programs) from R."""
        if self._r_fit_id is None:
            raise RuntimeError("EBMF is not fitted. Call prepare_anndata() first.")
        ro = _ensure_rpy2()
        from rpy2.robjects import numpy2ri, pandas2ri

        out = {}
        with (ro.default_converter + numpy2ri.converter + pandas2ri.converter).context():
            for i, ct in enumerate(self._scores):
                ref = f'{self._r_fit_id}_fit[["loadings"]][[{i + 1}]]'
                W = np.asarray(ro.r(ref))
                genes = [str(x) for x in ro.r(f"rownames({ref})")]
                cols = [str(x) for x in ro.r(f"colnames({ref})")]
                out[ct] = pd.DataFrame(W, index=genes, columns=cols)
        return out
