"""Fit MOFA-FLEX only (no patpy scoring); save Z per config.
Usage: <mofaflex-python> fit_mofaflex_only.py <arm_dir> <prefix> <zdir>
"""
import os, sys, time, warnings, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

CONFIGS = [{"weight_prior": "Normal", "n_factors": 10},
           {"weight_prior": "Horseshoe", "n_factors": 10}]

arm_dir, prefix, zdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(zdir, exist_ok=True)
if prefix == "mean_A":
    with open("/mnt/shared-workspace/combat_benchmark/pb_mean_A.pkl", "rb") as f:
        pb = pickle.load(f)
else:
    pb = {f[len(prefix)+2:-4]: pd.read_csv(os.path.join(arm_dir, f), index_col=0)
          for f in sorted(os.listdir(arm_dir)) if f.startswith(prefix + "__") and f.endswith(".csv")}
print(f"arm '{prefix}': {len(pb)} cell types", flush=True)

import torch, anndata as ad
import mofaflex as mf
torch.set_num_threads(8)

for cfg in CONFIGS:
    tag = f"mofaflex_{prefix}_{cfg['weight_prior']}_k{cfg['n_factors']}"
    if os.path.exists(f"{zdir}/{tag}_Z.csv"):
        print(f"skip {tag}", flush=True)
        continue
    t0 = time.time()
    views = {ct: ad.AnnData(df.astype(np.float32)) for ct, df in pb.items()}
    mo = mf.ModelOptions(n_factors=int(cfg["n_factors"]),
                         weight_prior=cfg["weight_prior"], init_factors="pca")
    to = mf.TrainingOptions(device="cpu", max_epochs=10000, seed=67,
                            early_stopper_patience=100, save_path=f"/workspace/mofaflex_h5/{tag}.h5")
    os.makedirs("/workspace/mofaflex_h5", exist_ok=True)
    model = mf.MOFAFLEX({"single_group": views}, mo, to)
    Z = model.get_factors(); Z = Z[list(Z.keys())[0]]
    Z.index = Z.index.astype(str)
    Z.to_csv(f"{zdir}/{tag}_Z.csv")
    print(f"{tag}: fit done ({round(time.time()-t0)}s), Z={Z.shape}", flush=True)
print("FIT DONE", flush=True)
