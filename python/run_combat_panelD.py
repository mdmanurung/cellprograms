"""Panel D (partial) for COMBAT: interpretability without pathway coherence.

- Loading sparsity: fraction of near-zero loadings (|L| < 0.01) and Gini
  concentration, per method (EBMF per cell type; MOFA per view).
- Factor-score technical associations: per factor, Spearman vs total cells
  per sample (library-size proxy) and kNN macro-F1 vs Institute/Pool
  (technical retention per factor).

Usage: python run_combat_panelD.py
"""
import os, sys, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
import combat_common as cc

warnings.filterwarnings("ignore", category=UserWarning)


def gini(x):
    x = np.abs(np.asarray(x, dtype=float))
    if x.sum() == 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def sparsity_stats(L):
    """L: genes x K. Returns (frac_near_zero, mean_gini_per_factor)."""
    if L is None or L.size == 0:
        return np.nan, np.nan
    frac = float((np.abs(L) < 0.01).mean())
    g = np.mean([gini(L[:, k]) for k in range(L.shape[1])])
    return frac, g


def factor_technical_associations(rep, st):
    """Per factor: Spearman vs cell counts; kNN macro-F1 vs institute/pool."""
    from patpy.tl import evaluate_representation
    D = ssd.squareform(ssd.pdist(rep.values, metric="euclidean"))
    common = rep.index.intersection(st.index)
    pos = rep.index.get_indexer(common)
    D = D[np.ix_(pos, pos)]
    st_al = st.loc[common]
    rows = []
    for f in rep.columns:
        row = {"factor": f}
        # per-factor distances: |score_i - score_j| (1-D representation)
        d1 = np.abs(rep.loc[common, f].values[:, None] - rep.loc[common, f].values[None, :])
        for tgt, task in [("institute", "classification"), ("pool", "classification")]:
            t = st_al[tgt]
            try:
                res = evaluate_representation(d1, t, method="knn", task=task)
                row[f"knn_{tgt}"] = res.get("score")
            except Exception:
                row[f"knn_{tgt}"] = np.nan
        if "total_cells" in st_al.columns:
            from scipy.stats import spearmanr
            v = rep.loc[common, f].values
            m = np.isfinite(st_al["total_cells"].values)
            if m.sum() > 10:
                row["spearman_total_cells"] = float(spearmanr(v[m], st_al["total_cells"].values[m]).statistic)
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    # total cells per sample (library-size proxy)
    adata, st2 = cc.load_combat()
    counts = adata.obs.groupby(cc.SAMPLE_KEY, observed=True).size()
    st["total_cells"] = counts.reindex(st.index)
    del adata

    # EBMF sparsity from full-fit loadings (panel C checkpoint)
    if cc.ckpt_exists("ebmf_full_fits.pkl"):
        full = cc.load_ckpt("ebmf_full_fits.pkl")
        rows = []
        for ct, info in full.items():
            L = info["loadings"]
            fnz, g = sparsity_stats(L)
            rows.append({"method": "EBMF", "view": ct, "K": 0 if L is None else L.shape[1],
                         "frac_near_zero": fnz, "mean_gini": g})
        pd.DataFrame(rows).to_csv("/mnt/results/combat_benchmark/panel_D_ebmf_sparsity.csv", index=False)
        print("EBMF sparsity saved", flush=True)

    # MOFA sparsity from hdf5 weights
    import h5py
    for tag in ["A", "B", "B2"]:
        p = f"{cc.WORK}/mofa_{tag}.hdf5"
        if not os.path.exists(p):
            continue
        rows = []
        with h5py.File(p, "r") as f:
            for view in f["expectations"]["W"].keys():
                W = f["expectations"]["W"][view][:]
                fnz, g = sparsity_stats(W)
                rows.append({"method": f"MOFAcellulaR-{tag}", "view": view,
                             "K": W.shape[1], "frac_near_zero": fnz, "mean_gini": g})
        pd.DataFrame(rows).to_csv(f"/mnt/results/combat_benchmark/panel_D_mofa{tag}_sparsity.csv", index=False)
        print(f"MOFA-{tag} sparsity saved", flush=True)

    # factor-level technical associations for each representation
    for name, f in [("EBMF", "rep_ebmf_A.pkl"), ("MOFAcellulaR", "rep_mofa_A.pkl"),
                    ("MOFAcellulaR-B", "rep_mofa_B.pkl"), ("MOFAcellulaR-B2", "rep_mofa_B2.pkl")]:
        if not cc.ckpt_exists(f):
            continue
        rep = cc.load_ckpt(f)
        fa = factor_technical_associations(rep, st)
        fa.insert(0, "method", name)
        fa.to_csv(f"/mnt/results/combat_benchmark/panel_D_{name.replace('-', '')}_factor_technical.csv", index=False)
        print(f"{name} factor-technical saved", flush=True)
