"""EBMF full-factorial tuning search on COMBAT (TMM+HVG arm).

32 configs = loading_prior(2) x max_factors(4) x var_type(2) x backfit(2).
Locked config (point_laplace/30/1/TRUE) is one of the 32. Per-config metrics:
disease kNN (corrected macro-F1) + WHO-ordinal Spearman. Checkpoint per config.

Usage: python run_tuning_search.py <arm_dir> <out_csv>
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import combat_common as cc

warnings.filterwarnings("ignore", category=UserWarning)

PRIORS = ["point_laplace", "point_normal"]
MAXF = [10, 15, 30, 50]
VART = [0, 1]
BACKFIT = [True, False]


def all_configs():
    out = []
    for p in PRIORS:
        for m in MAXF:
            for v in VART:
                for b in BACKFIT:
                    out.append({"loading_prior": p, "max_factors": m,
                                "var_type": v, "backfit": b})
    return out


def _worker(job):
    ct, df, cfg = job
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
    return ct, m.build_representation()


def fit_config(pb, cfg, workers=8):
    from concurrent.futures import ProcessPoolExecutor
    blocks = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for ct, block in ex.map(_worker, [(ct, pb[ct], cfg) for ct in sorted(pb)]):
            blocks.append(block)
    rep = pd.concat(blocks, axis=1).fillna(0.0)
    rep.columns = rep.columns.astype(str)
    return rep


def score_rep(rep, st):
    from patpy.tl import evaluate_representation
    import scipy.spatial.distance as ssd
    common = rep.index.intersection(st.index)
    pos = rep.index.get_indexer(common)
    D = ssd.squareform(ssd.pdist(rep.values[pos], metric="euclidean"))
    st_al = st.loc[common]
    out = {}
    t = st_al["disease"]
    mask = t.notna().values
    res = evaluate_representation(D[np.ix_(mask, mask)], t.dropna(), method="knn",
                                  task="classification")
    out["knn_disease"] = res["score"]
    t2 = st_al["who_ordinal"]
    mask2 = t2.notna().values
    res2 = evaluate_representation(D[np.ix_(mask2, mask2)], t2.dropna(), method="knn",
                                   task="regression")
    out["knn_who"] = res2["score"]
    out["n_factors"] = rep.shape[1]
    return out


def load_arm(arm_dir, prefix):
    pb = {}
    for f in sorted(os.listdir(arm_dir)):
        if f.startswith(prefix + "__") and f.endswith(".csv"):
            ct = f[len(prefix) + 2:-4]
            pb[ct] = pd.read_csv(os.path.join(arm_dir, f), index_col=0)
    return pb


if __name__ == "__main__":
    arm_dir, out_csv = sys.argv[1], sys.argv[2]
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    pb = load_arm(arm_dir, "tmm_hvg")
    print(f"arm loaded: {len(pb)} cell types", flush=True)

    done = set()
    if os.path.exists(out_csv):
        done = set(pd.read_csv(out_csv)[["loading_prior", "max_factors", "var_type",
                                         "backfit"]].astype(str).agg("|".join, axis=1))
    rows = []
    for i, cfg in enumerate(all_configs()):
        key = "|".join(str(v) for v in [cfg["loading_prior"], cfg["max_factors"],
                                        cfg["var_type"], cfg["backfit"]])
        tag = f"cfg{i:02d}_{cfg['loading_prior']}_f{cfg['max_factors']}_v{cfg['var_type']}_b{int(cfg['backfit'])}"
        if key in done:
            print(f"skip {tag} (done)", flush=True)
            continue
        t0 = time.time()
        try:
            rep = fit_config(pb, cfg)
            s = score_rep(rep, st)
            rows.append({**cfg, **s, "runtime_s": round(time.time() - t0), "config": tag})
            print(f"{tag}: disease={s['knn_disease']:.3f} who={s['knn_who']:.3f} "
                  f"K={s['n_factors']} ({s['runtime_s']}s)", flush=True)
        except Exception as e:
            rows.append({**cfg, "error": str(e)[:200], "config": tag})
            print(f"{tag}: ERROR {str(e)[:120]}", flush=True)
        # append-checkpoint (rewrite whole CSV; small)
        pd.DataFrame(rows).to_csv(out_csv, index=False) if not os.path.exists(out_csv) else \
            pd.concat([pd.read_csv(out_csv), pd.DataFrame(rows)]).to_csv(out_csv, index=False)
        rows = []
    print("SEARCH DONE", flush=True)
