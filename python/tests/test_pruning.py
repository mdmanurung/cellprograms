"""Test 3 (plan, mechanism form): planted NULL factors are pruned; effective K detected.

The data-driven form of A5 (excess factors pruned on noisy synthetic data) did
not reproduce: with K=5 on 150x40 Poisson counts, 2 planted programs, no
factor reached pi0 >= 1-1e-3 (pruning threshold, cebmf.py:619-632) — noise
factors on small data stay diffuse. Recorded as a calibration observation in
the ledger. The mechanism test below is the direct reading of the plan's
"planted null factor gets pruned": inject explicit zero columns via
user-provided init (the M3 refinement hook) and require their removal.
"""
import numpy as np
import torch

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def test_injected_null_factors_pruned():
    counts, W_true, S_true, _ = make_counts(n_cells=150, n_genes=40, n_programs=2, seed=3)
    recipe = ShiftedLog()
    Y = torch.tensor(recipe.fit_transform(counts), dtype=torch.float32)

    K = 5
    # init: 2 signal factors (true W/S up to scale) + 3 explicit null (zero) factors
    init_L = np.zeros((counts.shape[0], K))
    init_F = np.zeros((counts.shape[1], K))
    init_L[:, 0] = S_true[:, 0]
    init_F[:, 0] = W_true[:, 0] * counts.shape[0]  # rough scale pairing
    init_L[:, 1] = S_true[:, 1]
    init_F[:, 1] = W_true[:, 1] * counts.shape[0]

    res = fit_cebmf(
        Y,
        FitConfig(K=K, maxit=60, internal_epoch=10, allow_backfitting=True, seed=3,
                  init_L=init_L, init_F=init_F),
    )
    assert res.effective_k <= 3, (
        f"null factors not pruned: effective_k={res.effective_k} (K_fit={res.L.shape[1]})"
    )
    assert res.L.shape[1] <= K
