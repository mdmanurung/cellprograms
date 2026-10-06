"""MOFA-FLEX on Stephenson: sum pseudobulk -> TMM + scran HVG arm (winner arm).
Usage: <mofaflex-python> run_stephenson_mofaflex.py <h5ad> <zdir>
"""
import os, sys, subprocess, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, anndata as ad
warnings.filterwarnings("ignore")
MIN_CELLS = 10

h5ad_path, zdir = sys.argv[1], sys.argv[2]
sum_dir = "/workspace/stephenson_arms/sum"
os.makedirs(sum_dir, exist_ok=True); os.makedirs(zdir, exist_ok=True)

adata = ad.read_h5ad(h5ad_path)
sk, ck, dk = "sample_id", "author_cell_type", "Status"
keep = ~adata.obs[sk].isin(adata.obs.loc[adata.obs[dk] == "Non_covid", sk].unique())
adata = adata[keep.values].copy()
X = adata.X
if hasattr(X, "toarray"):
    X = X.toarray()
obs = adata.obs[[sk, ck]].copy(); obs[sk] = obs[sk].astype(str)
for ct, idx in obs.groupby(ck)[ck].indices.items():
    samp = obs[sk].values[idx]
    df = pd.DataFrame(X[idx], index=samp, columns=adata.var_names.astype(str))
    cnt = df.groupby(level=0).size()
    keep_s = cnt[cnt >= MIN_CELLS].index
    if len(keep_s) < 5:
        continue
    pb = df.loc[df.index.isin(keep_s)].groupby(level=0).sum()
    pb.to_csv(f"{sum_dir}/{ct}.csv")
print("sum pseudobulk written", flush=True)

r = subprocess.run(["Rscript", "/workspace/cellprograms/python/preprocess_arms.R",
                    sum_dir, sum_dir, "/workspace/stephenson_arms/processed_tmm"],
                   capture_output=True, text=True,
                   env={**os.environ, "R_HOME": os.popen("R RHOME").read().strip(),
                        "R_LIBS_USER": "/workspace/.Rlib"})
print(r.stdout[-300:], r.stderr[-300:], flush=True)
if r.returncode != 0:
    raise RuntimeError("preprocess failed")

import torch
torch.set_num_threads(8)
import mofaflex as mf
pb = {f[len("tmm_hvg__"):-4]: pd.read_csv(os.path.join("/workspace/stephenson_arms/processed_tmm", f), index_col=0)
      for f in sorted(os.listdir("/workspace/stephenson_arms/processed_tmm"))
      if f.startswith("tmm_hvg__") and f.endswith(".csv")}
print(f"arm: {len(pb)} cell types", flush=True)
for prior in ["Normal", "Horseshoe"]:
    tag = f"mofaflex_stephenson_tmmhvg_{prior}_k10"
    if os.path.exists(f"{zdir}/{tag}_Z.csv"):
        continue
    t0 = time.time()
    views = {ct.replace(".", "_"): ad.AnnData(df.astype(np.float32)) for ct, df in pb.items()}
    mo = mf.ModelOptions(n_factors=10, weight_prior=prior, init_factors="pca")
    to = mf.TrainingOptions(device="cpu", max_epochs=10000, seed=67,
                            early_stopper_patience=100, save_path=f"/workspace/mofaflex_h5/{tag}.h5")
    model = mf.MOFAFLEX({"single_group": views}, mo, to)
    Z = model.get_factors(); Z = Z[list(Z.keys())[0]]
    Z.index = Z.index.astype(str)
    Z.to_csv(f"{zdir}/{tag}_Z.csv")
    print(f"{tag}: fit done ({round(time.time()-t0)}s)", flush=True)
print("STEPHENSON MOFAFLEX DONE", flush=True)
