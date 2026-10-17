"""Workstream G / M2 benchmark: real-data projection validation + 100k memory profile.

(a) Real data (UCSF CD8 singlecell granularity, 29998 cells x 2000 HVGs):
    discovery fit on a 10k-cell subsample (dense), then chunked projection of
    ALL 29998 cells with F/tau frozen. Reports matched |cor| and relative
    shape error between the projection's L and the discovery fit's L on the
    shared 10k cells (projection-vs-joint discrepancy), plus wall-clock and
    peak memory for the 30k projection.
(b) Memory profile at 100k x 3000: synthetic sparse counts at ~89% zero
    (matching the real singlecell granularity), chunked projection with
    10k-row chunks; reports peak RSS. The 5M extrapolation uses this
    measurement (design doc 03 §5).
"""
import json
import resource
import time

import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.io import mmread
from scipy.stats import spearmanr

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog
from cebmf_ext.project import project_scores

ROOT = "/mnt/shared-workspace/shared/audit/F_granularity"
K = 10


def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def part_a():
    X = mmread(f"{ROOT}/gran/singlecell/counts.mtx").tocsr()
    meta = pd.read_csv(f"{ROOT}/gran/singlecell/obs_meta.tsv", sep="\t")
    genes = [line.rstrip("\n") for line in open(f"{ROOT}/gran/singlecell/genes.tsv")]
    # drop all-zero genes (66 in this granularity; same drop as the cNMF arm)
    nz = np.asarray((X > 0).sum(axis=0)).ravel() > 0
    X = X[:, nz]
    genes = [g for g, keep in zip(genes, nz) if keep]

    rng = np.random.default_rng(20261017)
    disc_idx = np.sort(rng.choice(X.shape[0], size=10000, replace=False))

    X_disc = X[disc_idx].toarray()
    recipe = ShiftedLog()
    Y_disc = torch.tensor(recipe.fit_transform(X_disc), dtype=torch.float32)
    t0 = time.time()
    res = fit_cebmf(Y_disc, FitConfig(K=K, maxit=100, internal_epoch=10, seed=20261017))
    t_disc = time.time() - t0
    F = res.F
    tau_vec = np.asarray(res.tau)
    tau = float(np.mean(1.0 / tau_vec)) if tau_vec.ndim == 2 else float(tau_vec)
    print(f"[disc] K_eff={res.effective_k} runtime={t_disc:.0f}s peak={peak_gb():.2f}GB")

    t0 = time.time()
    L_proj = project_scores(X, F, tau, n_sweeps=5, chunk_size=10000, transform=recipe)
    t_proj = time.time() - t0
    print(f"[proj] shape={L_proj.shape} runtime={t_proj:.0f}s peak={peak_gb():.2f}GB")

    # comparison on the discovery cells: joint-fit L vs projection L
    L_joint = res.L
    A = L_joint - L_joint.mean(0)
    B = L_proj[disc_idx] - L_proj[disc_idx].mean(0)
    C = np.abs(np.corrcoef(A.T, B.T)[: L_joint.shape[1], L_joint.shape[1]:])
    matched = C.max(axis=1)
    rels = []
    for k in range(L_joint.shape[1]):
        j = int(np.argmax(C[k]))
        if matched[k] < 0.5:
            rels.append(np.nan)
            continue
        a = A[:, k]
        b = L_proj[disc_idx][:, j] - L_proj[disc_idx][:, j].mean()
        sc = np.dot(a, b) / max(np.dot(b, b), 1e-300)
        rels.append(np.linalg.norm(a - sc * b) / max(np.linalg.norm(a), 1e-300))
    out = pd.DataFrame({
        "factor": range(L_joint.shape[1]),
        "match_cor_proj_vs_joint": matched,
        "rel_err_proj_vs_joint": rels,
    })
    out.to_csv(f"{ROOT}/G_m2_projection_validation.csv", index=False)
    print(out.round(4).to_string(index=False))

    summary = {
        "part": "a_real_data_projection",
        "n_cells_total": int(X.shape[0]),
        "n_genes_after_zero_drop": int(X.shape[1]),
        "n_discovery_cells": 10000,
        "K_eff_discovery": int(res.effective_k),
        "tau_scalar_used": tau,
        "t_discovery_s": t_disc,
        "t_projection_s": t_proj,
        "projection_chunk_size": 10000,
        "peak_mem_gb": peak_gb(),
        "median_match_cor": float(np.nanmedian(matched)),
        "median_rel_err": float(np.nanmedian(rels)),
    }
    return summary


def part_b(n_cells=100000, n_genes=3000, sparsity=0.89):
    rng = np.random.default_rng(1)
    nnz = int(n_cells * n_genes * (1 - sparsity))
    rows = rng.integers(0, n_cells, size=nnz)
    cols = rng.integers(0, n_genes, size=nnz)
    vals = rng.poisson(1.0, size=nnz).astype(np.float32)
    X = sparse.csr_matrix((vals, (rows, cols)), shape=(n_cells, n_genes))
    del rows, cols, vals

    recipe = ShiftedLog()
    # fit preprocessing state on a small dense sample (as in practice)
    sample = X[:5000].toarray()
    recipe.fit_transform(sample)

    F = np.abs(rng.normal(0, 0.5, size=(n_genes, K))).astype(np.float32)
    tau = 1.0
    t0 = time.time()
    L = project_scores(X, F, tau, n_sweeps=3, chunk_size=10000, transform=recipe)
    t_proj = time.time() - t0
    summary = {
        "part": "b_memory_profile_100k",
        "n_cells": n_cells, "n_genes": n_genes, "sparsity": sparsity,
        "chunk_size": 10000, "n_sweeps": 3,
        "t_projection_s": t_proj,
        "peak_mem_gb": peak_gb(),
        "L_shape": list(L.shape),
    }
    return summary


if __name__ == "__main__":
    s_a = part_a()
    s_b = part_b()
    with open(f"{ROOT}/G_m2_benchmark_summary.json", "w") as f:
        json.dump({"a": s_a, "b": s_b}, f, indent=1)
    print(json.dumps({"a": s_a, "b": s_b}, indent=1))
