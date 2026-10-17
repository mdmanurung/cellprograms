"""Test 1 (plan): determinism — two identical fits are bit-identical."""
import numpy as np
import torch

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def test_bit_identical_refits():
    counts, _, _, _ = make_counts(seed=1)
    Y = torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)
    cfg = FitConfig(K=6, maxit=30, internal_epoch=10, seed=1)
    r1 = fit_cebmf(Y, cfg)
    r2 = fit_cebmf(Y, cfg)
    assert torch.equal(torch.tensor(r1.L), torch.tensor(r2.L))
    assert torch.equal(torch.tensor(r1.F), torch.tensor(r2.F))
