"""Preprocessing attribution: top-5 EBMF configs x {lognorm_all, lognorm_hvg}.

Usage: python run_attribution.py <processed_dir> <out_csv>
"""
import os, sys, time, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import combat_common as cc
from run_tuning_search import all_configs, fit_config, score_rep, load_arm

warnings.filterwarnings("ignore")

TOP5 = [  # (prior, maxf, vart, backfit) from TMM+HVG search ranking
    ("point_laplace", 10, 0, False),
    ("point_normal", 10, 1, True),
    ("point_normal", 10, 0, True),
    ("point_laplace", 10, 0, True),
    ("point_laplace", 10, 1, True),
]

if __name__ == "__main__":
    processed_dir, out_csv = sys.argv[1], sys.argv[2]
    _, st = cc.load_combat()
    st.index = st.index.astype(str)
    arms = {"lognorm_all": cc.load_ckpt("pb_mean_A.pkl"),
            "lognorm_hvg": load_arm(processed_dir, "lognorm_hvg")}
    print({k: len(v) for k, v in arms.items()}, flush=True)

    rows, done = [], set()
    if os.path.exists(out_csv):
        d = pd.read_csv(out_csv)
        rows = d.to_dict("records")
        done = set(d.loc[d["knn_disease"].notna(), "config"])

    for arm, pb in arms.items():
        for j, (p, m, v, b) in enumerate(TOP5):
            cfg = {"loading_prior": p, "max_factors": m, "var_type": v, "backfit": b}
            tag = f"{arm}_{p}_f{m}_v{v}_b{int(b)}"
            if tag in done:
                print(f"skip {tag}", flush=True)
                continue
            t0 = time.time()
            try:
                rep = fit_config(pb, cfg)
                s = score_rep(rep, st)
                rows.append({**cfg, "arm": arm, **s,
                             "runtime_s": round(time.time() - t0), "config": tag})
                print(f"{tag}: disease={s['knn_disease']:.3f} who={s['knn_who']:.3f} "
                      f"K={s['n_factors']} ({round(time.time()-t0)}s)", flush=True)
            except Exception as e:
                import traceback; traceback.print_exc()
                rows.append({**cfg, "arm": arm, "error": str(e)[:200], "config": tag})
            pd.DataFrame(rows).to_csv(out_csv, index=False)  # rewrite whole file
    print("ATTRIBUTION DONE", flush=True)
