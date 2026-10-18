# Python helpers for cellprogramsMC, executed inside the basilisk
# environment via basilisk::basiliskRun(). Sourced with runpy so the
# functions below are callable from R.
#
# SEACells API (PyPI 0.3.3, verified against installed source):
#   model = SEACells.core.SEACells(ad, build_kernel_on='X_pca', n_SEACells=...,
#                                  use_gpu=False, n_waypoint_eigs=10,
#                                  n_neighbors=15, convergence_epsilon=1e-5,
#                                  use_sparse=True)
#   model.initialize_archetypes(); model.fit(max_iter=...)
#   model.get_hard_assignments()  -> DataFrame('SEACell'), index = cells
#   model.A_                      -> (n_metacells x n_cells) assignment matrix
#   SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell',
#                                      summarize_layer='raw') -> AnnData (mc x genes, summed counts)
#
# Preprocessing follows the hdWGCNA "Alternative metacell algorithms" vignette
# (which follows the SEACells tutorial): normalize_per_cell -> log1p -> HVG ->
# PCA, no scaling before PCA.

import os

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse


def run_seacells_group(mtx_path, cells_csv, genes_csv, out_dir,
                       n_seacells, n_top_genes=1500, n_pc=50,
                       pca_csv=None, n_waypoint_eigs=10, waypt_proportion=0.9,
                       n_neighbors=15, convergence_epsilon=1e-5,
                       max_iter=20, seed=20261016, use_sparse=True,
                       mode="sum"):
    """Run SEACells on one group and export metacell artifacts.

    Parameters mirror the R-side arguments. Exports (into ``out_dir``):
      metacell_mtx.mtx      genes x metacells, summed raw counts
      metacell_obs.csv      per-metacell metadata (n_cells, group size)
      soft_assignments.mtx  cells x metacells soft assignment (A_.T), optional
      fit_stats.csv         RSS trace summary and fit diagnostics
    """
    import scanpy as sc
    import SEACells

    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.RandomState(int(seed))  # R numerics arrive as float64
    np.random.seed(int(seed))

    # ---- load counts (genes x cells as written by R) -> cells x genes ----
    X_gxc = scipy.io.mmread(mtx_path)
    X_gxc = scipy.sparse.csr_matrix(X_gxc)
    cells = pd.read_csv(cells_csv)["cell"].astype(str).values
    genes = pd.read_csv(genes_csv)["gene"].astype(str).values
    if X_gxc.shape != (len(genes), len(cells)):
        raise ValueError(
            "counts matrix shape %s does not match genes/cells (%d, %d)"
            % (X_gxc.shape, len(genes), len(cells))
        )
    X = X_gxc.T.tocsr()
    # R-written mtx loads as integer dtype; in-place normalize_per_cell would
    # fail casting float64 output back into int64. Cast to float up front.
    X = X.astype(np.float32)

    ad = sc.AnnData(X)
    ad.obs_names = cells
    ad.var_names = genes
    # retain raw counts for aggregation
    raw_ad = sc.AnnData(X.copy())
    raw_ad.obs_names = ad.obs_names
    raw_ad.var_names = ad.var_names
    ad.raw = raw_ad

    # ---- preprocessing (SEACells tutorial / hdWGCNA vignette pattern) ----
    sc.pp.normalize_per_cell(ad)
    sc.pp.log1p(ad)
    n_hvg = int(min(n_top_genes, ad.n_vars - 1))
    sc.pp.highly_variable_genes(ad, n_top_genes=n_hvg)

    if pca_csv is not None:
        P = pd.read_csv(pca_csv, index_col=0)
        P.index = P.index.astype(str)
        P = P.loc[ad.obs_names]
        ad.obsm["X_pca"] = P.values.astype(np.float64)
    else:
        sc.tl.pca(ad, n_comps=int(min(n_pc, min(ad.shape) - 1)),
                  use_highly_variable=True)

    # ---- SEACells model ----
    n_seacells = int(min(n_seacells, ad.n_obs - 1))
    if n_seacells < 1:
        raise ValueError("group too small for SEACells (n_cells=%d)" % ad.n_obs)
    model = SEACells.core.SEACells(
        ad,
        build_kernel_on="X_pca",
        n_SEACells=n_seacells,
        use_gpu=False,
        verbose=True,
        n_waypoint_eigs=int(n_waypoint_eigs),
        n_neighbors=int(n_neighbors),
        convergence_epsilon=float(convergence_epsilon),
        # NOTE: SEACells 0.3.3 (PyPI GPU fork) has no use_sparse kwarg; the
        # parameter is kept in the signature for API compatibility but not
        # passed through (W5 scratch variant first documented this drift).
    )
    model.waypt_proportion = float(waypt_proportion)
    # SEACells 0.3.3 (GPU-accelerated fork) defers kernel construction:
    # __init__ sets K=None and initialize_archetypes/_fit require an explicit
    # construct_kernel_matrix() call first (upstream tutorial API drift).
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    converged = True
    try:
        model.fit(max_iter=int(max_iter))
    except RuntimeWarning as w:
        # this fork raises on non-convergence; record and continue — the RSS
        # trace in fit_stats documents how close the fit got
        converged = False
        print("SEACells fit did not converge: %s" % w)

    # ---- assignments ----
    hard = model.get_hard_assignments()  # DataFrame 'SEACell', index = cells
    ad.obs["SEACell"] = hard["SEACell"].values

    sead = SEACells.core.summarize_by_SEACell(
        ad, SEACells_label="SEACell", summarize_layer="raw"
    )  # metacells x genes, summed raw counts

    mc_sizes = ad.obs["SEACell"].value_counts()
    mc_names = list(sead.obs_names)
    sizes = np.array([int(mc_sizes.get(m, 0)) for m in mc_names])

    # metacell matrix: genes x metacells
    mc_gxm = scipy.sparse.csr_matrix(sead.X).T
    if mode == "average":
        mc_gxm = mc_gxm @ scipy.sparse.diags(1.0 / np.maximum(sizes, 1))
    scipy.io.mmwrite(os.path.join(out_dir, "metacell_mtx.mtx"), mc_gxm)

    obs = pd.DataFrame(
        {
            "n_cells": sizes,
            "n_metacells_requested": int(n_seacells),
            "n_cells_group": int(ad.n_obs),
        },
        index=mc_names,
    )
    obs.to_csv(os.path.join(out_dir, "metacell_obs.csv"))

    # soft assignments: A_ is (n_metacells x n_cells); export cells x metacells
    A = scipy.sparse.csc_matrix(model.A_.T)
    scipy.io.mmwrite(os.path.join(out_dir, "soft_assignments.mtx"), A)

    # fit statistics: RSS trace if available
    stats = {
        "n_cells": [int(ad.n_obs)],
        "n_genes": [int(ad.n_vars)],
        "n_hvg": [int(ad.var["highly_variable"].sum())],
        "n_pcs": [int(ad.obsm["X_pca"].shape[1])],
        "n_metacells": [int(len(mc_names))],
        "median_mc_size": [float(np.median(sizes))],
        "min_mc_size": [int(sizes.min())],
        "max_mc_size": [int(sizes.max())],
        "seed": [int(seed)],
    }
    rss = getattr(model, "RSS_iters", None)
    if rss is not None and len(rss) > 0:
        rss_arr = np.asarray(rss).reshape(-1)
        stats["rss_initial"] = [float(rss_arr[0])]
        stats["rss_final"] = [float(rss_arr[-1])]
        stats["rss_n_iter"] = [int(len(rss_arr))]
    stats["converged"] = [bool(converged)]
    pd.DataFrame(stats).to_csv(os.path.join(out_dir, "fit_stats.csv"), index=False)

    return {
        "n_metacells": int(len(mc_names)),
        "n_cells": int(ad.n_obs),
        "out_dir": out_dir,
    }
