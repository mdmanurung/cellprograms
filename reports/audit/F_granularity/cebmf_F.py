"""Workstream F (PILOT): cebmf_torch arms on three input granularities.

Same-model-family comparator to the R ebmf arm: identical shifted-log
preprocessing of COUNTS (y = log1p(x/(alpha*s_i)), alpha=1, library-size size
factors), then cEBMF with the closest prior/noise match to the package
operating point:
  prior_L = prior_F = 'laplace'  (point-Laplace ~ flashier point_laplace)
  noise_type = COLUMN_WISE       (~ var_type=1, per-gene variances)
  K = 10 initial, allow_backfitting = True (pi0-based factor pruning)
cebmf_torch treats NaN entries as missing (native). SVD init is deterministic,
so one run per granularity; a pseudobulk repeat verifies determinism.
"""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import torch
from scipy.io import mmread
from cebmf_torch import cEBMF
from cebmf_torch.cebmf.cebmf import NoiseType  # noqa: F401  (path pinned at clone rev)

ROOT = "/mnt/shared-workspace/shared/audit/F_granularity"
ARMS = ["pseudobulk", "metacell", "singlecell"]
K_INIT = 10
MAXIT = 100
SEED = 20261017


def shifted_log(counts):
    sf = counts.sum(axis=1)
    sf = sf / sf.mean()
    return np.log1p(counts / (1.0 * sf[:, None]))


def donor_scores(S, meta, arm):
    donor = meta["donor"].astype(str).values
    if arm == "pseudobulk":
        return pd.DataFrame(S, index=donor)
    if arm == "metacell":
        w = meta["n_cells"].values.astype(float)
        rows = []
        for d in pd.unique(donor):
            ii = donor == d
            rows.append(np.average(S[ii], axis=0, weights=w[ii]))
        return pd.DataFrame(rows, index=pd.unique(donor))
    rows = [S[donor == d].mean(axis=0) for d in pd.unique(donor)]
    return pd.DataFrame(rows, index=pd.unique(donor))


def run_arm(arm, torch_device):
    d = os.path.join(ROOT, "gran", arm)
    X = mmread(os.path.join(d, "counts.mtx")).toarray()
    meta = pd.read_csv(os.path.join(d, "obs_meta.tsv"), sep="\t")
    genes = [line.rstrip("\n") for line in open(os.path.join(d, "genes.tsv"))]

    Y = torch.tensor(shifted_log(X), dtype=torch.float32, device=torch_device)
    t0 = time.time()
    model = cEBMF(Y, K=K_INIT, prior_L="laplace", prior_F="laplace",
                  noise_type=NoiseType.COLUMN_WISE, allow_backfitting=True)
    fit = model.fit(maxit=MAXIT)
    el = time.time() - t0

    L = fit.L.cpu().numpy()  # obs x K
    F = fit.F.cpu().numpy()  # genes x K
    K_fit = L.shape[1]
    # surviving factors: non-degenerate columns after pruning
    col_normL = np.linalg.norm(L, axis=0)
    K_eff = int((col_normL > 1e-8).sum())
    obj = [float(x) for x in fit.history_obj] if fit.history_obj is not None else []

    W = pd.DataFrame(F, index=genes)
    S = pd.DataFrame(L, index=meta["id"].astype(str).values)
    W.to_csv(f"{ROOT}/F_cebmf_{arm}_loadings.csv")
    S.to_csv(f"{ROOT}/F_cebmf_{arm}_scores.csv")
    ds = donor_scores(L, meta, arm)
    ds.to_csv(f"{ROOT}/F_cebmf_{arm}_donor_scores.csv")
    pd.DataFrame([{"arm": arm, "K_init": K_INIT, "K_fit": K_fit, "K_eff": K_eff,
                   "maxit": MAXIT, "runtime_s": el, "seed": SEED,
                   "obj_first": obj[0] if obj else np.nan,
                   "obj_last": obj[-1] if obj else np.nan,
                   "obj_delta_last5": float(np.mean(np.diff(obj[-6:]))) if len(obj) > 5 else np.nan}]
                ).to_csv(f"{ROOT}/F_cebmf_{arm}_summary.csv", index=False)
    print(f"[{arm}] K_fit={K_fit} K_eff={K_eff} runtime={el:.0f}s obj_delta={obj[-1] - obj[-2] if len(obj) > 1 else 'NA'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="pseudobulk,metacell,singlecell")
    ap.add_argument("--cpu", action="store_true", default=True)
    args = ap.parse_args()
    dev = torch.device("cpu")
    torch.manual_seed(SEED)
    for arm in args.arms.split(","):
        run_arm(arm.strip(), dev)
    # determinism check on pseudobulk
    d = os.path.join(ROOT, "gran", "pseudobulk")
    X = mmread(os.path.join(d, "counts.mtx")).toarray()
    Y = torch.tensor(shifted_log(X), dtype=torch.float32, device=dev)
    fit2 = cEBMF(Y, K=K_INIT, prior_L="laplace", prior_F="laplace",
                 noise_type=NoiseType.COLUMN_WISE, allow_backfitting=True).fit(maxit=MAXIT)
    fit3 = cEBMF(Y, K=K_INIT, prior_L="laplace", prior_F="laplace",
                 noise_type=NoiseType.COLUMN_WISE, allow_backfitting=True).fit(maxit=MAXIT)
    same = bool(torch.allclose(fit2.L, fit3.L) and torch.allclose(fit2.F, fit3.F))
    pd.DataFrame([{"check": "cebmf_pseudobulk_determinism", "bit_identical": same}]).to_csv(
        f"{ROOT}/F_cebmf_determinism.csv", index=False)
    print("cebmf F arms done; determinism:", same)
