"""Workstream G / M1 API-verification pilot: covariate-moderated cEBMF (Approach B).

Question (user, 2026-10-17): celltype as covariate? Celltype labels are NOT
available in the built granularities (CD8-only cohort; builder kept donor id
only) — would require re-extraction from the source h5ad. Donor AGE is
available (gran/donor_age.tsv), external (B-L1 compliant), and on-mission
(immune-aging). This pilot runs age as the covariate; celltype pilot deferred
pending label re-extraction.

Design: pseudobulk granularity (99 donors x 2000 HVGs, shifted-log).
  Fit A (reference): laplace/laplace, covariate-free  (= workstream F arm).
  Fit B (moderated): prior_L="cash" with X_l = standardized age; prior_F laplace.
  Fit B2 (control):  prior_L="cash" with X_l = SHUFFLED age (permutation control).
Metrics: finiteness, factor matching A-vs-B (Hungarian |cor|), usage-vs-age
Spearman per factor (B vs A-post-hoc baseline: the same Spearman on A's usage
is the "post-hoc regression" comparator of design-doc criterion B3).
Status: API-verification pilot, NOT a registered moderation claim (B-L3
requires held-out donors; that is M5 scope).
"""
import json
import time

import numpy as np
import pandas as pd
import torch
from scipy.io import mmread
from scipy.stats import spearmanr
from cebmf_torch import cEBMF
from cebmf_torch.cebmf.cebmf import NoiseType

ROOT = "/mnt/shared-workspace/shared/audit/F_granularity"
K = 10
MAXIT = 100


def shifted_log(counts):
    sf = counts.sum(axis=1)
    sf = sf / sf.mean()
    return np.log1p(counts / (1.0 * sf[:, None]))


def main():
    X = mmread(f"{ROOT}/gran/pseudobulk/counts.mtx").toarray()
    meta = pd.read_csv(f"{ROOT}/gran/pseudobulk/obs_meta.tsv", sep="\t")
    age = pd.read_csv(f"{ROOT}/gran/donor_age.tsv", sep="\t").rename(
        columns={"donor_id": "donor"})
    df = meta.merge(age, on="donor", how="left")
    assert df["age"].notna().all(), "missing ages"
    Y = torch.tensor(shifted_log(X), dtype=torch.float32)
    age_z = torch.tensor(
        ((df["age"].values - df["age"].values.mean()) / df["age"].values.std()),
        dtype=torch.float32,
    ).reshape(-1, 1)

    rng = np.random.default_rng(20261017)
    age_shuf = age_z[torch.tensor(rng.permutation(len(age_z)))]

    fits = {}
    for name, kwargs in [
        ("A_laplace", dict(prior_L="laplace", prior_F="laplace")),
        ("B_cash_age", dict(prior_L="cash", prior_F="laplace", X_l=age_z)),
        ("B2_cash_shuffled", dict(prior_L="cash", prior_F="laplace", X_l=age_shuf)),
    ]:
        t0 = time.time()
        m = cEBMF(Y, K=K, noise_type=NoiseType.COLUMN_WISE,
                  allow_backfitting=True, **kwargs)
        fit = m.fit(maxit=MAXIT)
        L = fit.L.cpu().numpy()
        obj = [float(x) for x in fit.history_obj]
        fits[name] = L
        print(f"[{name}] K_fit={L.shape[1]} finite={np.isfinite(L).all()} "
              f"obj_finite={np.isfinite(obj).all()} runtime={time.time()-t0:.0f}s")

    rows = []
    for name, L in fits.items():
        for k in range(L.shape[1]):
            r = spearmanr(df["age"].values, L[:, k]).statistic
            rows.append({"fit": name, "factor": k, "spearman_usage_age": r})
    out = pd.DataFrame(rows)
    out.to_csv(f"{ROOT}/G_pilot_cash_age_usage_spearman.csv", index=False)

    # cross-fit matching A vs B
    def best_match_cors(A, B):
        Ac = A - A.mean(0); Bc = B - B.mean(0)
        C = np.corrcoef(Ac.T, Bc.T)[:A.shape[1], A.shape[1]:]
        C = np.abs(C)
        vals = []
        Bc_used = set()
        for i in np.argsort(-C.max(axis=1)):
            j = int(np.argmax(C[i]))
            if j in Bc_used:
                order = np.argsort(-C[i])
                j = next((int(x) for x in order if x not in Bc_used), -1)
            if j >= 0:
                Bc_used.add(j); vals.append(float(C[i, j]))
        return vals

    mc = best_match_cors(fits["A_laplace"], fits["B_cash_age"])
    pd.DataFrame({"n_matched": [len(mc)], "match_cors": [";".join(f"{v:.3f}" for v in mc)]}).to_csv(
        f"{ROOT}/G_pilot_cash_age_match.csv", index=False)
    print("A-vs-B matched |cor|:", [round(v, 3) for v in mc])

    summary = {
        "pilot": "cash_age_moderation_api_verification",
        "K_init": K, "maxit": MAXIT,
        "n_donors": int(len(df)),
        "age_available": True,
        "celltype_available": False,
        "celltype_note": "granularities carry donor id only; CD8-only cohort; "
                         "subtype labels need source-h5ad re-extraction",
        "fits": {n: {"K_fit": int(L.shape[1]), "finite": bool(np.isfinite(L).all())}
                 for n, L in fits.items()},
    }
    with open(f"{ROOT}/G_pilot_cash_age_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print("pilot done")


if __name__ == "__main__":
    main()
