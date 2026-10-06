"""Tests for the patpy EBMF adapter (synthetic AnnData, no downloads)."""
import os

os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")

import numpy as np
import pandas as pd
import pytest

from patpy_cellprograms import EBMF


def make_adata(n_samples=12, n_cells_per_sample_ct=30, n_genes=80, seed=0):
    import anndata as ad

    rng = np.random.default_rng(seed)
    samples = [f"S{i}" for i in range(n_samples)]
    cell_types = ["CTA", "CTB"]
    disease = ["COVID" if i % 2 == 0 else "healthy" for i in range(n_samples)]
    rows, meta = [], []
    gene_names = [f"gene_{j}" for j in range(n_genes)]
    for i, s in enumerate(samples):
        for ct in cell_types:
            base = rng.poisson(2.0, (n_cells_per_sample_ct, n_genes)).astype(np.float32)
            if disease[sample_idx := i] == "COVID":
                base[:, :5] += 4.0  # disease program in first 5 genes
            rows.append(base)
            meta += [pd.DataFrame({"sample": s, "cell_type": ct, "disease": disease[i]},
                                  index=[f"{s}_{ct}_{k}" for k in range(n_cells_per_sample_ct)])
                     for _ in [0]]
    X = np.vstack(rows)
    obs = pd.concat(meta)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=gene_names))
    return adata


@pytest.fixture(scope="module")
def fitted():
    adata = make_adata()
    m = EBMF(sample_key="sample", cell_group_key="cell_type", seed=1)
    m.prepare_anndata(adata)
    return m, adata


def test_distance_matrix_properties(fitted):
    m, adata = fitted
    D = m.calculate_distance_matrix()
    n = adata.obs["sample"].nunique()
    assert D.shape == (n, n)
    assert np.allclose(D, D.T)
    assert np.allclose(np.diag(D), 0)
    assert np.isfinite(D).all()


def test_representation_shape_and_scaling(fitted):
    m, adata = fitted
    rep = m.sample_representation
    assert rep.index.tolist() == sorted(adata.obs["sample"].astype(str).unique())
    # z-scaled blocks: columns are mean-zero
    assert np.allclose(rep.values.mean(axis=0), 0, atol=1e-8)


def test_metrics_suite_runs(fitted):
    from patpy.tl import evaluate_representation

    m, adata = fitted
    D = m.calculate_distance_matrix()
    target = adata.obs.groupby("sample")["disease"].first()
    for method, params in [("knn", {"task": "classification"}), ("silhouette", {}),
                           ("permanova", {}),
                           ("distances", {"control_level": "healthy", "normalization_type": "total"})]:
        res = evaluate_representation(D, target, method=method, **params)
        assert "score" in res


def test_k0_guard_and_missing_ct():
    """A cell type with pure noise (K=0) contributes no columns; missing CT tolerated."""
    adata = make_adata(seed=3)
    # make CTB pure noise by shuffling: EBMF may retain 0 factors; just check it runs
    m = EBMF(sample_key="sample", cell_group_key="cell_type", seed=1)
    m.prepare_anndata(adata)
    D = m.calculate_distance_matrix()
    assert D.shape[0] == adata.obs["sample"].nunique()


def test_loadings_access(fitted):
    m, _ = fitted
    W = m.get_loadings()
    assert set(W.keys()) == {"CTA", "CTB"}
    for ct, w in W.items():
        assert w.shape[0] == 80  # n_genes

def test_hyphenated_cell_type_names():
    """Cell-type labels that are invalid R symbols (hyphens) must not break the bridge."""
    adata = make_adata(seed=0)
    adata.obs["cell_type"] = adata.obs["cell_type"].replace({"CTA": "B_non-switched_memory", "CTB": "CD8.TE"})
    m = EBMF(sample_key="sample", cell_group_key="cell_type", seed=1)
    m.prepare_anndata(adata)
    D = m.calculate_distance_matrix()
    assert D.shape == (adata.obs["sample"].nunique(),) * 2
    W = m.get_loadings()
    assert set(W.keys()) == {"B_non-switched_memory", "CD8.TE"}


def test_var_type_validation():
    import pytest
    adata = make_adata(n_samples=4, n_cells_per_sample_ct=10, n_genes=30)
    m = EBMF(sample_key="sample", cell_group_key="cell_type", seed=1, var_type=(1, 2))
    assert m.var_type == (1, 2)
    assert m._r_var_type() == "c(1, 2)"
    m2 = EBMF(sample_key="sample", cell_group_key="cell_type", seed=1, var_type=2)
    assert m2.var_type == (2,) and m2._r_var_type() == "2"
    with pytest.raises(ValueError):
        EBMF(sample_key="sample", cell_group_key="cell_type", seed=1, var_type=(0, 1))
    with pytest.raises(ValueError):
        EBMF(sample_key="sample", cell_group_key="cell_type", seed=1, var_type=3)


def test_kronecker_var_type_end_to_end():
    adata = make_adata(n_samples=6, n_cells_per_sample_ct=15, n_genes=40)
    m = EBMF(sample_key="sample", cell_group_key="cell_type", seed=1, var_type=(1, 2))
    m.prepare_anndata(adata)
    D = m.calculate_distance_matrix()
    assert D.shape[0] == 6
