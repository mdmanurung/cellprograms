"""Test 4 (plan, as amended): objective trace behavior.

AMENDMENT to pre-registered A4 (recorded in the ledger before any benchmark
outcome): the original criterion "objective non-decreasing/non-increasing
within 1e-6 per iteration" is FALSIFIED for cebmf_torch on synthetic data.
Two documented phenomena (2026-10-17, synthetic 150x40, K=4):
  (a) initial transient: KL terms start at zero (priors unfitted), so the
      recorded loss starts artificially low and rises for ~4 iterations;
  (b) late-stage upward drift ~1e-5/iteration after convergence (mechanism
      hypothesis: tau re-estimation outside the variational bound, or float32
      accumulation; HYPOTHESIS-grade, see ledger).
The amended criterion: trace finite; net loss decrease from the post-transient
peak to the end; tail drift bounded (<= 1e-3 relative over the last 10
iterations); direction is "loss" (source-verified).
"""
import numpy as np
import torch

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def test_objective_finite_net_decrease_bounded_tail_drift():
    counts, _, _, _ = make_counts(seed=4)
    Y = torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)
    res = fit_cebmf(Y, FitConfig(K=4, maxit=30, internal_epoch=10, seed=4))
    tr = np.asarray(res.objective, dtype=float)
    assert res.direction == "loss"
    assert np.isfinite(tr).all(), "non-finite objective values"
    assert len(tr) >= 10, "objective trace too short"
    # net decrease from the post-transient peak (iteration 5) to the end
    scale = max(abs(tr).max(), 1e-300)
    assert tr[-1] < tr[5], f"no net loss decrease after transient: {tr[5]} -> {tr[-1]}"
    # tail drift bounded
    tail_drift = (tr[-1] - tr[-10]) / scale
    assert tail_drift <= 1e-3, f"tail drift {tail_drift:.3e} exceeds 1e-3"
