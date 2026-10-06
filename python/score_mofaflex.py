"""Score saved MOFA-FLEX Z matrices. Usage: python score_mofaflex.py <zdir> <out_csv>"""
import os, sys, warnings
os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import combat_common as cc
from run_tuning_search import score_rep
warnings.filterwarnings("ignore")

zdir, out_csv = sys.argv[1], sys.argv[2]
_, st = cc.load_combat(); st.index = st.index.astype(str)
rows = []
for f in sorted(os.listdir(zdir)):
    if f.endswith("_Z.csv"):
        tag = f[:-6]
        Z = pd.read_csv(os.path.join(zdir, f), index_col=0)
        Z.index = Z.index.astype(str)
        s = score_rep(Z, st)
        parts = tag.split("_")  # mofaflex <prefix> <prior> k<K>
        rows.append({"weight_prior": parts[2], "n_factors": parts[3][1:], **s,
                     "runtime_s": "", "config": tag})
        print(f"{tag}: disease={s['knn_disease']:.3f} who={s['knn_who']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(out_csv, index=False)
print("SCORE DONE", flush=True)
