"""Result contract: everything downstream (R cellprograms included) consumes.

Design rule (audit 01 §2): fields the backend cannot provide are serialized as
null — never manufactured. cebmf_torch exposes point estimates only, so
lfsr/posterior sd are null by construction.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import torch

SCHEMA_VERSION = "cebmf-ext-contract/0.1"

REQUIRED_NULLABLE = ["lfsr_loadings", "lfsr_scores", "posterior_sd_loadings", "posterior_sd_scores"]


def save_result(
    out_dir: str,
    name: str,
    result,  # cebmf_ext.fit.FitResult
    gene_ids: list[str],
    obs_ids: list[str],
    preprocessing: dict | None = None,
) -> dict:
    """Write loadings.csv (genes x K), scores.csv (obs x K), contract.json."""
    os.makedirs(out_dir, exist_ok=True)
    K = result.L.shape[1]
    cols = [f"{name}_{k}" for k in range(K)]
    pd.DataFrame(result.F, index=gene_ids, columns=cols).to_csv(
        os.path.join(out_dir, f"{name}_loadings.csv"))
    pd.DataFrame(result.L, index=obs_ids, columns=cols).to_csv(
        os.path.join(out_dir, f"{name}_scores.csv"))

    import cebmf_torch

    contract = {
        "schema_version": SCHEMA_VERSION,
        "name": name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "n_obs": int(result.L.shape[0]),
        "n_genes": int(result.F.shape[0]),
        "K_fit": int(K),
        "effective_k": int(result.effective_k),
        "converged": bool(result.converged),
        "objective_direction": result.direction,
        "objective_trace": result.objective,
        "runtime_s": result.runtime_s,
        "peak_mem_gb": result.peak_mem_gb,
        "config": result.config.to_dict(),
        "preprocessing": preprocessing or {},
        "backend": {
            "package": "cebmf_torch",
            "torch_version": torch.__version__,
            "cebmf_torch_path": os.path.dirname(cebmf_torch.__file__),
            "device": "cpu",
        },
        # Capability gaps of the backend, explicit and null (audit 01 §1):
        "lfsr_loadings": None,
        "lfsr_scores": None,
        "posterior_sd_loadings": None,
        "posterior_sd_scores": None,
    }
    with open(os.path.join(out_dir, f"{name}_contract.json"), "w") as f:
        json.dump(contract, f, indent=1)
    return contract


def validate_contract(contract: dict) -> list[str]:
    """Return a list of schema violations (empty = valid)."""
    errs = []
    for key in ["schema_version", "name", "n_obs", "n_genes", "K_fit", "effective_k",
                "converged", "objective_trace", "runtime_s", "config", "backend"]:
        if key not in contract:
            errs.append(f"missing required key: {key}")
    for key in REQUIRED_NULLABLE:
        if key not in contract:
            errs.append(f"missing nullable capability key (must be present, may be null): {key}")
    if "preprocessing" not in contract:
        errs.append("missing preprocessing provenance")
    if contract.get("K_fit") is not None and contract.get("effective_k") is not None:
        if contract["effective_k"] > contract["K_fit"]:
            errs.append("effective_k exceeds K_fit")
    return errs


def load_result(out_dir: str, name: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    with open(os.path.join(out_dir, f"{name}_contract.json")) as f:
        contract = json.load(f)
    W = pd.read_csv(os.path.join(out_dir, f"{name}_loadings.csv"), index_col=0)
    S = pd.read_csv(os.path.join(out_dir, f"{name}_scores.csv"), index_col=0)
    return contract, W, S
