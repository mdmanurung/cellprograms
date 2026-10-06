"""Reconstruction + sparsity for the TUNED winner config (pl/10/v1/backfit).
Fixes run_winner_characterize.py, which used ebmf_full_fit defaults (K=30).
Usage: python run_winner_recon.py <processed_dir> <out_csv>
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import combat_common as cc
from run_tuning_search import load_arm
from run_combat_panelC import reconstruction_metrics

warnings.filterwarnings("ignore")
CFG = {"loading_prior": "point_laplace", "max_factors": 10, "var_type": 1, "backfit": True}


def _fit_one(item):
    ct, df, cfg = item
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
    os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from patpy_cellprograms import EBMF
    m = EBMF(sample_key="sample", cell_group_key="ct", seed=67,
             loading_prior=cfg["loading_prior"], max_factors=int(cfg["max_factors"]),
             var_type=int(cfg["var_type"]), backfit=bool(cfg["backfit"]))
    m._fit_in_r({ct: df})
    m._fitted = True
    block = m.build_representation()
    Ldf = m.get_loadings().get(ct)
    return ct, block, (None if Ldf is None else (Ldf.values, Ldf))


def gini(x):
    x = np.abs(np.asarray(x, dtype=float)).ravel()
    x = x[x > 0]
    if x.size == 0:
        return np.nan
    x = np.sort(x)
    n = x.size
    return float((2 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())) - (n + 1) / n)


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    processed_dir, out_csv = sys.argv[1], sys.argv[2]
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    arms = {"lognorm_all": cc.load_ckpt("pb_mean_A.pkl"),
            "lognorm_hvg": load_arm(processed_dir, "lognorm_hvg"),
            "tmm_hvg": load_arm(processed_dir, "tmm_hvg")}
    rows = []
    for arm, pb in arms.items():
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=8) as ex:
            for ct, block, L in ex.map(_fit_one, [(ct, pb[ct], CFG) for ct in sorted(pb)]):
                Lvals = None if L is None else L[0]
                Y = pb[ct].values.astype(np.float64)
                if Lvals is not None and Lvals.shape[0] != Y.shape[1]:
                    Lvals = L[1].reindex(pb[ct].columns, fill_value=0.0).values
                S = block.reindex(pb[ct].index).fillna(0.0).values
                r2, r2h = reconstruction_metrics(Y, Lvals, S)
                rows.append({"arm": arm, "cell_type": ct,
                             "K": 0 if Lvals is None else Lvals.shape[1],
                             "r2": r2, "r2_heldout": r2h,
                             "gini": np.nan if Lvals is None else gini(Lvals)})
        print(f"{arm}: recon done ({time.time()-t0:.0f}s)", flush=True)
    rec = pd.DataFrame(rows)
    rec.to_csv(out_csv, index=False)
    print(rec.groupby("arm")[["K", "r2", "r2_heldout", "gini"]].mean().round(3).to_string(), flush=True)
    print("RECON DONE", flush=True)
