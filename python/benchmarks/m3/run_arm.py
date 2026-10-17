"""M3 arm runner (pre-registered 04-experiment-plan.md §1, as amended by 04a).

One process = one (arm, rate, rep) cell. All arms share the registered shifted-log
recipe fit on TRAIN-donor cells only; all arms fit on the 80% training donors;
held-out prediction projects the 20% test donors with F fixed (donor-safe, B-L2).

Arms: S1 from-scratch | S2 pseudobulk-init+project | S3 metacell-init+project |
S4pb/S4mc init+refine (from S2/S3 factors) | S5 cNMF anchor (run_cnmf.py).
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import torch

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.project import project_scores
from m3_metrics import fixed_entry_mask, heldout_prediction, hungarian_match, match_cor_matrix, top_jaccard
from m3_norm import make_recipe
from m3_synth import generate

SEED0 = 20261016
MAXIT_FULL, MAXIT_REFINE = 100, 50
N_METACELL_PER_DONOR = 5


def fit_with_fallback(Y, cfg):
    """Fit with recorded K-halving fallback on backend divergence.

    The confirmed silent-divergence regime (small matrix + over-fitted K; ledger,
    commit 5452201) makes K=10 discovery fits on pseudobulk/metacell matrices
    (40-200 rows) diverge. The adapter raises loudly; we retry at K/2 and record
    the fallback - this quantifies the divergence's operational cost.
    """
    try:
        res = fit_cebmf(Y, cfg)
        return res, False
    except RuntimeError:
        cfg2 = FitConfig(**{**cfg.to_dict(), "K": max(cfg.K // 2, 2), "init_L": None, "init_F": None})
        res = fit_cebmf(Y, cfg2)
        return res, True


def tau_scalar(res) -> float:
    tau = np.asarray(res.tau)
    return float(np.mean(1.0 / tau)) if tau.ndim == 2 else float(tau)


def project(counts_csr, F, tau, recipe, chunk=2500, sweeps=5):
    return project_scores(counts_csr, F, tau, n_sweeps=sweeps, chunk_size=chunk, transform=recipe)


def aggregate_donor_means(Y: np.ndarray, donor_ids: np.ndarray) -> np.ndarray:
    """Mean normalized expression per TRAINING donor (donor ids present in Y)."""
    return np.stack([Y[donor_ids == d].mean(0) for d in np.unique(donor_ids)])


def metacell_matrix(Y: np.ndarray, donor_ids: np.ndarray,
                    k_per_donor: int, seed: int):
    from sklearn.cluster import KMeans
    rows, weights = [], []
    for d in np.unique(donor_ids):
        idx = np.where(donor_ids == d)[0]
        km = KMeans(n_clusters=k_per_donor, n_init=5, random_state=seed).fit(Y[idx])
        for c in range(k_per_donor):
            m = idx[km.labels_ == c]
            rows.append(Y[m].mean(0))
            weights.append(len(m))
    M = np.asarray(rows)
    w = np.asarray(weights, dtype=np.float64)[:, None]
    return M, w  # metacell x genes, weights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True,
                    choices=["S1", "S2", "S3", "S4pb", "S4mc"])
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--rep", type=int, required=True)
    ap.add_argument("--beta", type=float, default=0.0)
    ap.add_argument("--norm", default="shiftedlog", choices=["shiftedlog", "pflog", "clr"],
                    help="Amendment A4: normalization variant (S1/S3 sensitivity arms)")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    t_start = time.time()
    g = generate(args.rep, args.rate, args.beta)
    counts = g["counts"]
    train_idx = np.where(~g["test_mask"])[0]
    test_idx = np.where(g["test_mask"])[0]
    X_train = counts[train_idx].toarray()
    X_test = counts[test_idx].toarray()

    recipe = make_recipe(args.norm)
    Y_train = recipe.fit_transform(X_train)  # recipe fit on TRAIN donors only
    Y_test_true = recipe.transform(X_test)
    Yt = torch.tensor(Y_train, dtype=torch.float32)
    seed = SEED0 + args.rep
    zero_rate = float((X_train == 0).mean())

    row = dict(arm=args.arm, rate=args.rate, rep=args.rep, beta=args.beta,
               n_train=len(train_idx), n_test=len(test_idx), zero_rate=zero_rate,
               seed=seed, norm=args.norm)
    if hasattr(recipe, "alpha_hat") and recipe.alpha_hat is not None:
        row["alpha_hat"] = round(recipe.alpha_hat, 6)

    if args.arm in ("S2", "S3", "S4pb", "S4mc"):
        if args.arm in ("S2", "S4pb"):
            agg = aggregate_donor_means(Y_train, g["donor_ids"][train_idx])
            init_arm = "S2"
        else:
            agg, w = metacell_matrix(Y_train, g["donor_ids"][train_idx],
                                     N_METACELL_PER_DONOR, seed)
            init_arm = "S3"
        t0 = time.time()
        res_init, div_init = fit_with_fallback(torch.tensor(agg, dtype=torch.float32),
                                               FitConfig(K=10, maxit=MAXIT_FULL, internal_epoch=10, seed=seed))
        t_init = time.time() - t0
        F_init = res_init.F
        L_proj = project(counts, F_init, tau_scalar(res_init), recipe)
        if args.arm in ("S4pb", "S4mc"):
            t0 = time.time()
            res, div = fit_with_fallback(Yt, FitConfig(K=10, maxit=MAXIT_REFINE, internal_epoch=10,
                                                       seed=seed, init_L=L_proj[train_idx],
                                                       init_F=F_init))
            t_fit = time.time() - t0
            L_train = res.L
            L_test = project(counts[test_idx], res.F, tau_scalar(res), recipe)
            row.update(init_arm=init_arm, t_init_s=t_init, diverged_init=div_init, diverged=div)
        else:
            res = res_init
            t_fit = t_init
            L_train, L_test = L_proj[train_idx], L_proj[test_idx]
            row.update(diverged_init=div_init, diverged=False)
    else:  # S1
        t0 = time.time()
        res, div = fit_with_fallback(Yt, FitConfig(K=10, maxit=MAXIT_FULL, internal_epoch=10, seed=seed))
        t_fit = time.time() - t0
        L_train = res.L
        L_test = project(counts[test_idx], res.F, tau_scalar(res), recipe)
        row["diverged"] = div

    matches = hungarian_match(L_train, g["Z_true"][train_idx])
    marker_sets = [set(m) for m in g["markers"]]
    for m in matches:
        t = m["truth_col"]
        m["top50_jaccard"] = top_jaccard(res.F, m["est_col"], marker_sets[t])
    # Per-rep linear ceiling: best top-10 SVD |cor| per truth column on the SAME
    # normalized train matrix (Amendment A3: gates are evaluated relative to this
    # ceiling because compositional coupling caps absolute recovery for ANY linear
    # method; absolute bars are also reported at aggregation time).
    U, S, _ = np.linalg.svd(Y_train - Y_train.mean(0), full_matrices=False)
    C_svd = match_cor_matrix(U[:, :10] * S[:10], g["Z_true"][train_idx])
    svd_ceiling = C_svd.max(axis=0).tolist()
    mask = fixed_entry_mask(Y_test_true.shape, args.rep)
    row.update(heldout_prediction(L_test, res.F, Y_test_true, mask))
    row.update(
        k_eff=int(res.effective_k), converged=bool(res.converged),
        t_fit_s=round(t_fit, 2), t_total_s=round(time.time() - t_start, 2),
        peak_mem_gb=round(res.peak_mem_gb, 3),
        matches=matches, svd_ceiling=[round(c, 4) for c in svd_ceiling],
        final_objective=float(res.objective[-1]) if len(res.objective) else None,
    )
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"{args.arm}_rate{args.rate}_rep{args.rep:02d}.json")
    with open(out, "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: v for k, v in row.items() if k != "matches"}, indent=None))
    print("matches:", [(m["truth_col"], round(m["match_cor"], 3),
                         round(m["top50_jaccard"], 2)) for m in matches])


if __name__ == "__main__":
    main()
