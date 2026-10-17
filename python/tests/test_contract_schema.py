"""Test 5 (plan): contract schema — save/validate roundtrip; capability gaps are null."""
import numpy as np
import torch

from cebmf_ext import FitConfig, fit_cebmf, save_result, validate_contract
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def test_contract_roundtrip(tmp_path):
    counts, _, _, _ = make_counts(seed=1)  # seed=5 data diverges at K=4 (see test_divergence.py)
    recipe = ShiftedLog()
    Y = torch.tensor(recipe.fit_transform(counts), dtype=torch.float32)
    res = fit_cebmf(Y, FitConfig(K=4, maxit=20, internal_epoch=10, seed=1))  # seed=5/K=4 diverges (see test_divergence.py)
    contract = save_result(
        str(tmp_path), "test", res,
        gene_ids=[f"GENE{i}" for i in range(counts.shape[1])],
        obs_ids=[f"CELL{i}" for i in range(counts.shape[0])],
        preprocessing=recipe.to_dict(),
    )
    errs = validate_contract(contract)
    assert errs == [], f"contract violations: {errs}"
    for key in ["lfsr_loadings", "lfsr_scores", "posterior_sd_loadings", "posterior_sd_scores"]:
        assert contract[key] is None, f"{key} must be null (backend provides no uncertainty)"
    assert contract["effective_k"] <= contract["K_fit"]
    assert contract["preprocessing"]["recipe"] == "shifted_log"
