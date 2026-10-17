# Workstream F (PILOT): cNMF arms on three input granularities.
# Reuses the workstream-E cNMF pattern (counts in, own pipeline, 100 restarts,
# own density-silhouette K selection). Pilot deviations from E, documented in
# TASKS.md: K grid widened to 2..15 (real tissue, unknown K; E was sims with
# K_true=3); num_highvar_genes=None (input genes are already the top-2000 HVGs
# computed in build_granularity.py); single fixed seed 20261017 per arm
# (stability comes from cNMF's internal restarts + the ebmf seed trio).
import argparse
import os
import shutil
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from cnmf import cNMF
from scipy.io import mmread

warnings.filterwarnings("ignore")


def cnmf_silhouette(obj, k):
    st = obj.consensus(k=k, skip_density_and_return_after_stats=True,
                       show_clustering=False, close_clustergram_fig=True).stats
    v = st.iloc[0] if isinstance(st, pd.DataFrame) else st
    return float(v["silhouette"]), float(v["prediction_error"])


def run_arm(arm, gran_root, out_root, kmin, kmax, n_restarts, seed, work_root=None):
    d = os.path.join(gran_root, arm)
    X = mmread(os.path.join(d, "counts.mtx")).tocsr()
    meta = pd.read_csv(os.path.join(d, "obs_meta.tsv"), sep="\t")
    genes = [line.rstrip("\n") for line in open(os.path.join(d, "genes.tsv"))]
    counts = pd.DataFrame(X.toarray(), index=meta["id"].astype(str).values, columns=genes)
    # Drop all-zero gene columns (cNMF normalization divides by gene totals -> 0/0 NaN;
    # singlecell subsample zeroed 66/2000 HVGs). Recorded in the summary row.
    nz = counts.columns[(counts.values > 0).any(axis=0)]
    n_zero_genes = counts.shape[1] - len(nz)
    if n_zero_genes:
        counts = counts[nz]

    # cNMF writes h5ad/npz intermediates with random access — must be on local
    # disk (/workspace), NOT the S3-backed shared mount.
    work = os.path.join(work_root or "/workspace/cnmf_F_work", f"cnmf_work_{arm}")
    os.makedirs(work, exist_ok=True)
    h5 = os.path.join(work, f"{arm}.h5ad")
    ad.AnnData(X=counts.values.astype(np.float32),
               obs=pd.DataFrame(index=counts.index),
               var=pd.DataFrame(index=counts.columns)).write_h5ad(h5)
    obj = cNMF(output_dir=work, name=arm)
    components = list(range(kmin, kmax + 1))
    t0 = time.time()
    obj.prepare(counts_fn=h5, components=components, seed=seed, n_iter=n_restarts,
                num_highvar_genes=counts.shape[1])  # all input genes (already top-2000 HVGs); None breaks HVG selection
    obj.factorize()
    obj.combine(components=components)
    t_fit = time.time() - t0

    sil = {}
    for k in components:
        try:
            sil[k] = cnmf_silhouette(obj, k)
        except Exception as e:
            print(f"[{arm}] silhouette failed at K={k}: {e}")
            sil[k] = (float("nan"), float("nan"))
    valid = {k: v for k, v in sil.items() if np.isfinite(v[0])}
    if not valid:
        raise RuntimeError(f"[{arm}] no valid silhouette across K grid")
    k_hat = max(valid, key=lambda k: (valid[k][0], -k))

    obj.consensus(k=k_hat, show_clustering=False, close_clustergram_fig=True)
    usage, spectra, _, _ = obj.load_results(K=k_hat, density_threshold=0.5)

    usage.to_csv(os.path.join(out_root, f"F_cnmf_{arm}_usage.csv"))
    spectra.to_csv(os.path.join(out_root, f"F_cnmf_{arm}_spectra.csv"))
    pd.DataFrame([{"arm": arm, "k_hat": k_hat, "n_obs": counts.shape[0], "n_genes": counts.shape[1],
                   "n_zero_genes_dropped": n_zero_genes,
                   "t_fit_s": t_fit, "seed": seed, "n_restarts": n_restarts,
                   "silhouette_khat": valid[k_hat][0],
                   "silhouette_by_k": ";".join(f"{k}:{valid[k][0]:.4f}" for k in sorted(valid))}]
                ).to_csv(os.path.join(out_root, f"F_cnmf_{arm}_summary.csv"), index=False)
    print(f"[{arm}] k_hat={k_hat} t_fit={t_fit:.0f}s sil={valid[k_hat][0]:.4f}")
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="pseudobulk,metacell,singlecell")
    ap.add_argument("--gran-root", default="/mnt/shared-workspace/shared/audit/F_granularity/gran")
    ap.add_argument("--out-root", default="/mnt/shared-workspace/shared/audit/F_granularity")
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=15)
    ap.add_argument("--restarts", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--work-root", default="/workspace/cnmf_F_work")
    args = ap.parse_args()
    for arm in args.arms.split(","):
        run_arm(arm.strip(), args.gran_root, args.out_root, args.kmin, args.kmax,
                args.restarts, args.seed, args.work_root)
