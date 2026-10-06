"""Cross-dataset validation: tuned EBMF winner (pl/10/v1/backfit) on Stephenson.

Builds log-norm mean pseudobulk + scran top-2000 HVGs per cell type
(same recipe as the COMBAT lognorm_hvg arm), fits the winner, scores
disease kNN + severity Spearman, bootstrap CI.

Usage: python run_stephenson_winner.py <h5ad> <outdir>
"""
import os, sys, subprocess, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import anndata as ad
import scipy.spatial.distance as ssd
from run_patpy_comparison import SEVERITY_ORDER, bootstrap_ci
from run_tuning_search import fit_config

warnings.filterwarnings("ignore")
CFG = {"loading_prior": "point_laplace", "max_factors": 10, "var_type": 1, "backfit": True}
MIN_CELLS = 10

if __name__ == "__main__":
    h5ad_path, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    mean_dir = "/workspace/stephenson_arms/mean"
    os.makedirs(mean_dir, exist_ok=True)

    adata = ad.read_h5ad(h5ad_path)
    sk, ck, dk = "sample_id", "author_cell_type", "Status"
    keep = ~adata.obs[sk].isin(adata.obs.loc[adata.obs[dk] == "Non_covid", sk].unique())
    adata = adata[keep.values].copy()
    print(f"cells: {adata.shape}, samples: {adata.obs[sk].nunique()}", flush=True)

    per_sample = adata.obs.groupby(sk).agg(
        status=(dk, lambda s: s.dropna().iloc[0]),
        severity=("Worst_Clinical_Status",
                  lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan))
    per_sample.index = per_sample.index.astype(str)

    # mean pseudobulk of log-normalized X per cell type
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    obs = adata.obs[[sk, ck]].copy()
    obs[sk] = obs[sk].astype(str)
    n_written = 0
    for ct, idx in obs.groupby(ck)[ck].indices.items():
        samp = obs[sk].values[idx]
        df = pd.DataFrame(X[idx], index=samp, columns=adata.var_names.astype(str))
        cnt = df.groupby(level=0).agg(["size"]).iloc[:, 0] if False else df.groupby(level=0).size()
        keep_s = cnt[cnt >= MIN_CELLS].index
        if len(keep_s) < 5:
            continue
        pb = df.loc[df.index.isin(keep_s)].groupby(level=0).mean()
        pb.to_csv(f"{mean_dir}/{ct}.csv")
        n_written += 1
    print(f"pseudobulk written: {n_written} cell types", flush=True)

    # scran HVG selection (same as preprocess_arms.R lognorm path)
    r = subprocess.run(["Rscript", "/workspace/cellprograms/python/preprocess_arms.R",
                        mean_dir, mean_dir, "/workspace/stephenson_arms/processed"],
                       capture_output=True, text=True,
                       env={**os.environ, "R_HOME": os.popen("R RHOME").read().strip(),
                            "R_LIBS_USER": "/workspace/.Rlib"})
    print(r.stdout[-500:], r.stderr[-500:], flush=True)
    if r.returncode != 0:
        raise RuntimeError("preprocess failed")

    from run_tuning_search import load_arm, score_rep
    pb = load_arm("/workspace/stephenson_arms/processed", "lognorm_hvg")
    print(f"arm loaded: {len(pb)} cell types", flush=True)
    t0 = time.time()
    rep = fit_config(pb, CFG)
    rep.index = rep.index.astype(str)
    print(f"fit done ({time.time()-t0:.0f}s), rep={rep.shape}", flush=True)

    common = rep.index.intersection(per_sample.index)
    pos = rep.index.get_indexer(common)
    D = ssd.squareform(ssd.pdist(rep.values[pos], metric="euclidean"))
    st_al = per_sample.loc[common]
    rows = []
    t = st_al["status"]
    mask = t.notna().values
    from patpy.tl import evaluate_representation
    res = evaluate_representation(D[np.ix_(mask, mask)], t.dropna(), method="knn", task="classification")
    lo, hi, mean_b = bootstrap_ci(D[np.ix_(mask, mask)], t.dropna(), "knn", {"task": "classification"})
    rows.append({"metric": "knn_disease", "score": res["score"], "boot_mean": mean_b,
                 "ci_low": lo, "ci_high": hi})
    print(f"knn_disease: {res['score']:.3f} boot {mean_b:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
    sev = st_al["severity"].map({v: i for i, v in enumerate(SEVERITY_ORDER)})
    mask2 = sev.notna().values
    res2 = evaluate_representation(D[np.ix_(mask2, mask2)], sev.dropna(), method="knn", task="regression")
    rows.append({"metric": "knn_severity", "score": res2["score"]})
    print(f"knn_severity: {res2['score']:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(f"{outdir}/stephenson_winner.csv", index=False)
    rep.to_csv(f"{outdir}/stephenson_winner_rep.csv")
    print("STEPHENSON DONE", flush=True)
