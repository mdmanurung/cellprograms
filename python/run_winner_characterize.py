"""Winner characterization on COMBAT: tuned EBMF (pl/10/v1/backfit) x 3 arms.

Full biology panel (A+B) with bootstrap CIs, reconstruction R2/held-out,
sparsity (Gini of |loadings|). Saves reps + per-arm summary CSV.

Usage: python run_winner_characterize.py <processed_dir> <outdir>
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
import combat_common as cc
from run_tuning_search import fit_config, load_arm
from run_combat_metrics import PANEL_TARGETS
from run_patpy_comparison import run_metrics, bootstrap_ci
from run_combat_panelC import reconstruction_metrics, ebmf_full_fit

warnings.filterwarnings("ignore")

CFG = {"loading_prior": "point_laplace", "max_factors": 10, "var_type": 1, "backfit": True}
N_BOOT = 200


def gini(x):
    x = np.abs(np.asarray(x, dtype=float)).ravel()
    x = x[x > 0]
    if x.size == 0:
        return np.nan
    x = np.sort(x)
    n = x.size
    return float((2 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())) - (n + 1) / n)


def dist_from_rep(rep, st):
    common = rep.index.intersection(st.index)
    pos = rep.index.get_indexer(common)
    return ssd.squareform(ssd.pdist(rep.values[pos], metric="euclidean")), st.loc[common]


if __name__ == "__main__":
    processed_dir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    arms = {"lognorm_all": cc.load_ckpt("pb_mean_A.pkl"),
            "lognorm_hvg": load_arm(processed_dir, "lognorm_hvg"),
            "tmm_hvg": load_arm(processed_dir, "tmm_hvg")}

    panel_rows, recon_rows = [], []
    for arm, pb in arms.items():
        t0 = time.time()
        rep = fit_config(pb, CFG)
        rep.index = rep.index.astype(str)
        cc.save_ckpt(rep, f"rep_tuned_{arm}.pkl")
        D, st_al = dist_from_rep(rep, st)
        targets = PANEL_TARGETS["A"](st_al)
        panel_rows += run_metrics(D, targets, f"EBMF-tuned-{arm}")
        targets_b = PANEL_TARGETS["B"](st_al)
        panel_rows += run_metrics(D, targets_b, f"EBMF-tuned-{arm}")
        for metric, task in [("knn_disease", "classification"), ("knn_who_ordinal", "regression")]:
            t = targets[metric][0]
            lo, hi, mean = bootstrap_ci(D, t, "knn", {"task": task}, n_boot=N_BOOT)
            panel_rows.append({"method": f"EBMF-tuned-{arm}", "metric": f"{metric}_boot",
                               "eval_method": "bootstrap", "score": mean,
                               "ci_low": lo, "ci_high": hi})
        print(f"{arm}: panels done ({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(panel_rows).to_csv(f"{outdir}/winner_panel.csv", index=False)

    for arm, pb in arms.items():
        full = ebmf_full_fit(pb, workers=8)
        for ct, info in full.items():
            L = info["loadings"]
            Y = pb[ct].values.astype(np.float64)
            if L is not None and L.shape[0] != Y.shape[1]:
                Ldf = info.get("loadings_df")
                L = Ldf.reindex(pb[ct].columns, fill_value=0.0).values
            S = info["block"].reindex(pb[ct].index).fillna(0.0).values
            r2, r2h = reconstruction_metrics(Y, L, S)
            recon_rows.append({"arm": arm, "cell_type": ct,
                               "K": 0 if L is None else L.shape[1],
                               "r2": r2, "r2_heldout": r2h,
                               "gini": np.nan if L is None else gini(L)})
        print(f"{arm}: recon done", flush=True)
    rec = pd.DataFrame(recon_rows)
    rec.to_csv(f"{outdir}/winner_reconstruction.csv", index=False)
    print(rec.groupby("arm")[["r2", "r2_heldout", "gini"]].mean().round(3).to_string(), flush=True)
    print("CHARACTERIZATION DONE", flush=True)
