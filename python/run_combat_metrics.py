"""COMBAT metrics: Panels A (biology) and B (technical retention).

Reuses the proven patpy evaluate_representation + bootstrap machinery from
run_patpy_comparison.py. Distance matrices come from representation checkpoints.

Usage: python run_combat_metrics.py <A|B> <regime>
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
import combat_common as cc
from run_patpy_comparison import run_metrics, bootstrap_ci

warnings.filterwarnings("ignore", category=UserWarning)

PANEL_TARGETS = {
    "A": lambda st: {
        "knn_disease": (st["disease"], "knn", {"task": "classification"}),
        "knn_source": (st["source"], "knn", {"task": "classification"}),
        "knn_who_ordinal": (st["who_ordinal"], "knn", {"task": "regression"}),
        "knn_outcome": (st["outcome"], "knn", {"task": "regression"}),
        "silhouette_disease": (st["disease"], "silhouette", {}),
        "permanova_disease": (st["disease"], "permanova", {"permutations": 999}),
        "distances_disease": (st["disease"], "distances",
                              {"control_level": "Healthy", "normalization_type": "total"}),
    },
    "B": lambda st: {
        "knn_institute": (st["institute"], "knn", {"task": "classification"}),
        "knn_pool": (st["pool"], "knn", {"task": "classification"}),
        "silhouette_pool": (st["pool"], "silhouette", {}),
    },
}

REPS = {
    "EBMF": "rep_ebmf_{regime}.pkl",
    "MOFAcellulaR": "rep_mofa_{regime}.pkl",
    "MOFAcellulaR-B2": "rep_mofa_{regime}2.pkl",
    "PCA": "rep_pca_{regime}.pkl",
    "Pseudobulk": "rep_pseudobulk_{regime}.pkl",
    "GroupedPseudobulk": "rep_gpb_{regime}.pkl",
    "RandomVector": "rep_rand_{regime}.pkl",
}


def main(panel, regime):
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    all_rows = []
    for name, pat in REPS.items():
        f = pat.format(regime=regime)
        if not cc.ckpt_exists(f):
            continue
        obj = cc.load_ckpt(f)
        if isinstance(obj, dict):  # precomputed distances (GroupedPseudobulk)
            rep = pd.DataFrame(obj["distances"], index=obj["samples"], columns=obj["samples"])
        else:
            rep = obj
        if rep.values.ndim == 3:
            raise ValueError(f"{name}: 3-D representation not supported")
        if rep.index.equals(rep.columns):  # already a distance matrix
            D = rep.values.astype(float)
            common = rep.index.intersection(st.index)
            pos = rep.index.get_indexer(common)
            D = D[np.ix_(pos, pos)]
        else:
            common = rep.index.intersection(st.index)
            pos = rep.index.get_indexer(common)
            D = ssd.squareform(ssd.pdist(rep.values[pos], metric="euclidean"))
        st_al = st.loc[common]
        targets = PANEL_TARGETS[panel](st_al)
        print(f"== {name} ({rep.shape[0]} samples x {rep.shape[1]} features) ==", flush=True)
        all_rows += run_metrics(D, targets, name)
        # bootstrap CIs for the headline targets only
        for tname in (["knn_disease", "knn_source", "knn_who_ordinal"] if panel == "A"
                      else ["knn_institute", "knn_pool"]):
            if tname not in targets:
                continue
            target, method, params = targets[tname]
            mask = target.notna().values
            tv = target.dropna()
            pos2 = np.where(mask)[0]
            lo, hi, mean_b = bootstrap_ci(D[np.ix_(pos2, pos2)], tv.values, method, params)
            all_rows.append({"method": name, "metric": f"{tname}_boot", "eval_method": method,
                             "score": mean_b, "ci_low": lo, "ci_high": hi,
                             "n_observations": int(mask.sum())})
            print(f"  {tname} boot: {mean_b:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
    out = pd.DataFrame(all_rows)
    os.makedirs("/mnt/results/combat_benchmark", exist_ok=True)
    path = f"/mnt/results/combat_benchmark/panel_{panel}_regime{regime}.csv"
    out.to_csv(path, index=False)
    print("saved", path, flush=True)
    print(out.to_string())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "A",
         sys.argv[2] if len(sys.argv) > 2 else "A")
