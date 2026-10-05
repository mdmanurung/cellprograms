"""Run the PATpy metric suite comparison on Stephenson PBMC:
EBMF (cellprograms) vs MOFA vs Pseudobulk vs CellGroupComposition vs RandomVector.

Usage: python run_patpy_comparison.py <adata.h5ad> <outdir>
"""
from __future__ import annotations

import json
import os
import sys
import time

os.environ.setdefault("R_HOME", os.popen("R RHOME").read().strip())
os.environ.setdefault("R_LIBS_USER", "/workspace/.Rlib")

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from patpy_cellprograms import EBMF

from patpy.tl import evaluate_representation
from patpy.tl.sample_representation import (
    CellGroupComposition,
    MOFA,
    Pseudobulk,
    RandomVector,
)

RNG = np.random.default_rng(2026)
N_BOOTSTRAP = 200
BOOTSTRAP_FRACTION = 0.8

SEVERITY_ORDER = ["Healthy", "Asymptomatic", "Mild", "Moderate", "Severe", "Critical ", "Death"]


def pseudobulk_target(adata, sample_key, col):
    """Per-sample value of an obs column (first non-NA per sample)."""
    return adata.obs.groupby(sample_key)[col].agg(lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan)


def run_metrics(distances, targets, label):
    rows = []
    for target_name, (target, method, params) in targets.items():
        t0 = time.time()
        try:
            res = evaluate_representation(distances, target, method=method, **params)
            rows.append({
                "method": label, "metric": target_name, "eval_method": method,
                "score": res.get("score"), "p_value": res.get("p_value"),
                "n_observations": res.get("n_observations"), "runtime_s": round(time.time() - t0, 1),
            })
            print(f"  {target_name}: {res.get('score')}", flush=True)
        except Exception as e:  # noqa: BLE001
            rows.append({"method": label, "metric": target_name, "eval_method": method,
                         "score": np.nan, "error": str(e)[:200], "runtime_s": round(time.time() - t0, 1)})
            print(f"  {target_name}: ERROR {str(e)[:120]}", flush=True)
    return rows


def bootstrap_ci(distances, target, method, params, n_boot=N_BOOTSTRAP, frac=BOOTSTRAP_FRACTION):
    n = distances.shape[0]
    k = max(10, int(round(frac * n)))
    scores = []
    for _ in range(n_boot):
        idx = RNG.choice(n, size=k, replace=False)
        sub = distances[np.ix_(idx, idx)]
        t_sub = np.asarray(target)[idx]
        try:
            res = evaluate_representation(sub, t_sub, method=method, **params)
            scores.append(res["score"])
        except Exception:  # noqa: BLE001
            continue
    if not scores:
        return np.nan, np.nan, np.nan
    return float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5)), float(np.mean(scores))


def main(adata_path, outdir):
    import anndata as ad

    os.makedirs(outdir, exist_ok=True)
    adata = ad.read_h5ad(adata_path)
    print("adata:", adata.shape, flush=True)

    # Stephenson schema; drop Non-covid (other infections) for a clean Covid-vs-Healthy task.
    sample_key, cell_group_key, disease_key = "sample_id", "author_cell_type", "Status"
    keep = ~adata.obs[sample_key].isin(
        adata.obs.loc[adata.obs[disease_key] == "Non_covid", sample_key].unique()
    )
    adata = adata[keep.values].copy()
    print("after dropping Non_covid samples:", adata.shape, "| samples:", adata.obs[sample_key].nunique(), flush=True)

    per_sample = adata.obs.groupby(sample_key).agg(
        status=(disease_key, lambda s: s.dropna().iloc[0]),
        severity=("Worst_Clinical_Status", lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan),
    )
    targets_full = {
        "knn_disease": (per_sample["status"], "knn", {"task": "classification"}),
        "silhouette_disease": (per_sample["status"], "silhouette", {}),
        "permanova_disease": (per_sample["status"], "permanova", {"permutations": 999}),
        "distances_disease": (per_sample["status"], "distances",
                              {"control_level": "Healthy", "normalization_type": "total"}),
        "knn_severity": (per_sample["severity"].map({v: i for i, v in enumerate(SEVERITY_ORDER)}),
                         "knn", {"task": "regression"}),
    }

    methods = {
        "EBMF": lambda: EBMF(sample_key=sample_key, cell_group_key=cell_group_key, seed=67),
        "MOFA": lambda: MOFA(sample_key=sample_key, cell_group_key=cell_group_key, seed=67, n_factors=10),
        "Pseudobulk": lambda: Pseudobulk(sample_key=sample_key, cell_group_key=cell_group_key, seed=67),
        "CellGroupComposition": lambda: CellGroupComposition(sample_key=sample_key, cell_group_key=cell_group_key, seed=67),
        "RandomVector": lambda: RandomVector(sample_key=sample_key, cell_group_key=cell_group_key, seed=67),
    }

    all_rows, boot_rows = [], []
    for name, factory in methods.items():
        print(f"\n=== {name} ===", flush=True)
        t0 = time.time()
        try:
            m = factory()
            m.prepare_anndata(adata)
            D = m.calculate_distance_matrix()
            print(f"fit+distances: {time.time() - t0:.0f}s, D={D.shape}", flush=True)
            all_rows.extend(run_metrics(D, targets_full, name))
            # Bootstrap CIs for the primary metric (knn disease).
            target = targets_full["knn_disease"][0]
            mask = target.notna().values
            lo, hi, mean_b = bootstrap_ci(D[np.ix_(mask, mask)], target.dropna().values, "knn",
                                          {"task": "classification"})
            boot_rows.append({"method": name, "metric": "knn_disease", "mean": mean_b,
                              "ci_low": lo, "ci_high": hi, "n_boot": N_BOOTSTRAP})
            print(f"{name} knn_disease bootstrap = {mean_b:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
        except Exception as e:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            all_rows.append({"method": name, "metric": "FAILED", "error": str(e)[:300]})
        pd.DataFrame(all_rows).to_csv(os.path.join(outdir, "metrics_full.csv"), index=False)
        pd.DataFrame(boot_rows).to_csv(os.path.join(outdir, "metrics_bootstrap.csv"), index=False)

    print("\nDone.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
