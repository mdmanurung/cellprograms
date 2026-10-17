"""Regression test: backend divergence must be detected LOUDLY.

Confirmed backend finding (ledger 2026-10-17): on 150x40 synthetic counts
(seed=2 data, K=6, laplace/laplace, COLUMN_WISE), the cebmf_torch objective
diverges — loss climbs for 4 iterations then tau becomes NaN at iteration 5 —
silently (no warning, no error; NaN factors returned). The adapter must raise.
"""
import numpy as np
import pytest
import torch

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def test_divergence_raises():
    counts, _, _, _ = make_counts(n_cells=150, n_genes=40, n_programs=3, seed=2)
    Y = torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)
    with pytest.raises(RuntimeError, match="non-finite"):
        fit_cebmf(Y, FitConfig(K=6, maxit=40, internal_epoch=10, seed=2))
