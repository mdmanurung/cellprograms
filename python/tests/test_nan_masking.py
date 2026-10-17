"""Test 2 (plan): NaN masking — masked fit statistically recovers the complete fit.

Config note: seed=2 with K=6 DIVERGES in the backend (see test_divergence.py);
the masking equivalence uses a non-diverging config. The divergence itself is
a separate, confirmed backend finding (ledger 2026-10-17).
"""
import numpy as np
import torch

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def _match_cor(S1, S2):
    """Best finite |cor| per non-degenerate column of S1 against columns of S2."""
    K1, K2 = S1.shape[1], S2.shape[1]
    out = np.zeros(K1)
    for i in range(K1):
        a = S1[:, i]
        if np.std(a) < 1e-12:
            out[i] = 0.0
            continue
        best = 0.0
        for j in range(K2):
            b = S2[:, j]
            if np.std(b) < 1e-12:
                continue
            r = np.corrcoef(a, b)[0, 1]
            if np.isfinite(r):
                best = max(best, abs(r))
        out[i] = best
    return out


def test_nan_masking_equivalence():
    counts, _, S_true, _ = make_counts(n_cells=150, n_genes=40, n_programs=3, seed=1)
    recipe = ShiftedLog()
    Y_full = recipe.fit_transform(counts)

    rng = np.random.default_rng(2)
    Y_mask = Y_full.copy()
    drop = rng.random(Y_mask.shape) < 0.10
    Y_mask[drop] = np.nan

    cfg = FitConfig(K=4, maxit=40, internal_epoch=10, seed=1)
    r_full = fit_cebmf(torch.tensor(Y_full, dtype=torch.float32), cfg)
    r_mask = fit_cebmf(torch.tensor(Y_mask, dtype=torch.float32), cfg)

    cors = _match_cor(r_full.L, r_mask.L)
    # Criterion (A2 form): the PLANTED programs (3) must be well-matched between
    # the complete and masked fits. Excess factors (K=4 > 3 planted) are allowed
    # to differ — they fit noise and are not required to be stable.
    n_well_matched = int((cors >= 0.8).sum())
    assert n_well_matched >= 3, (
        f"only {n_well_matched}/3 planted programs well-matched under 10% masking: {cors}"
    )
