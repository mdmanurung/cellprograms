"""Run the PATpy metric suite comparison: EBMF (cellprograms) vs built-in methods.

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


def pseudobulk_target(adata, sample_key, col):
    """Per-sample value of an obs column (first non-NA per sample)."""
    return adata.obs.groupby(sample_key)[col].agg(lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan)


def run_metrics(distances, targets, label):
    """Run the metric suite on one distance matrix against several targets."""
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
        except Exception as e:  # noqa: BLE001 - record and continue
            rows.append({"method": label, "metric": target_name, "eval_method": method,
                         "score": np.nan, "error": str(e)[:200], "runtime_s": round(time.time() - t0, 1)})
    return rows


def bootstrap_ci(distances, target, method, params, n_boot=N_BOOTSTRAP, frac=BOOTSTRAP_FRACTION):
    """Bootstrap score distribution over donor subsets (resample without replacement)."""
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
    print("adata:", adata.shape)

    # Column detection (Stephenson processed schema).
    sample_key = next(c for c in adata.obs.columns if c.lower() in ("sample", "patient", "donor", "patient_id", "sample_id"))
    cell_group_key = next(c for c in adata.obs.columns if c.lower() in ("cell_type", "celltype", "annotation", "cell_types"))
    disease_key = next(c for c in adata.obs.columns if c.lower() in ("disease", "status", "outcome", "condition"))
    print(f"sample_key={sample_key} cell_group_key={cell_group_key} disease_key={disease_key}")

    targets_full = {
        "knn_disease": (pseudobulk_target(adata, sample_key, disease_key), "knn", {"task": "classification"}),
        "silhouette_disease": (pseudobulk_target(adata, sample_key, disease_key), "silhouette", {}),
        "permanova_disease": (pseudobulk_target(adata, sample_key, disease_key), "permanova", {"permutations": 999}),
        "distances_disease": (pseudobulk_target(adata, sample_key, disease_key), "distances",
                              {"control_level": sorted(pseudobulk_target(adata, sample_key, disease_key).dropna().unique())[0],
                               "normalization_type": "total"}),
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
            rows = run_metrics(D, targets_full, name)
            all_rows.extend(rows)
            # Bootstrap CIs for the primary metric (knn disease).
            target = targets_full["knn_disease"][0]
            mask = target.notna().values
            lo, hi, mean_b = bootstrap_ci(D[np.ix_(mask, mask)], target.dropna().values, "knn",
                                          {"task": "classification"})
            boot_rows.append({"method": name, "metric": "knn_disease", "mean": mean_b,
                              "ci_low": lo, "ci_high": hi, "n_boot": N_BOOTSTRAP})
            print(f"{name} knn_disease = {mean_b:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"{name} FAILED: {e}", flush=True)
            all_rows.append({"method": name, "metric": "FAILED", "error": str(e)[:300]})
        pd.DataFrame(all_rows).to_csv(os.path.join(outdir, "metrics_full.csv"), index=False)
        pd.DataFrame(boot_rows).to_csv(os.path.join(outdir, "metrics_bootstrap.csv"), index=False)

    print("\nDone.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
