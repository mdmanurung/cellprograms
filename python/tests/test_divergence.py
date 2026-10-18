"""Regression test: backend divergence must be detected LOUDLY.

Confirmed backend finding (ledger 2026-10-17): on 150x40 synthetic counts
(seed=2 data, K=6, laplace/laplace, COLUMN_WISE), the cebmf_torch objective
diverges — loss climbs then tau becomes NaN — silently (no warning, no error;
NaN factors returned).

R1 contract (ledger 2026-10-18): the adapter retries a registered seed family
(SEED_FAMILY_ATTEMPTS=3, stride 1000) and re-raises with the attempt
diagnostics only if every family member diverges. A successful fallback emits
a RuntimeWarning and records the attempt log in
FitResult.guard_diagnostics["seed_family"].
"""
import warnings

import numpy as np
import pytest
import torch

import cebmf_ext.fit as fit_mod
from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog

from conftest import make_counts


def _diverging_input():
    counts, _, _, _ = make_counts(n_cells=150, n_genes=40, n_programs=3, seed=2)
    return torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)


def test_divergence_raises_when_family_disabled(monkeypatch):
    """With the seed family disabled (single attempt), divergence raises."""
    Y = _diverging_input()
    monkeypatch.setattr(fit_mod, "SEED_FAMILY_ATTEMPTS", 1)
    with pytest.raises(RuntimeError, match="non-finite"):
        fit_cebmf(Y, FitConfig(K=6, maxit=40, internal_epoch=10, seed=2))


def test_divergence_family_exhaustion_raises():
    """Registered diverging fixture (ledger 2026-10-17): every seed-family
    member diverges on this input, so the adapter re-raises with the attempt
    diagnostics attached (no silent NaN factors)."""
    Y = _diverging_input()
    with pytest.raises(RuntimeError, match="diverged for all 3 seed-family attempts"):
        fit_cebmf(Y, FitConfig(K=6, maxit=40, internal_epoch=10, seed=2))


def test_seed_family_fallback_warns_and_records(monkeypatch):
    """Fallback contract: a diverging first attempt (simulated) is retried with
    the next family seed; recovery emits a RuntimeWarning and records the log."""
    counts, _, _, _ = make_counts(n_cells=80, n_genes=60, n_programs=2, seed=1)
    Y = torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)
    real = fit_mod._fit_once
    state = {"n": 0}

    def flaky(Y, cfg):
        state["n"] += 1
        if state["n"] == 1:
            raise RuntimeError("cEBMF objective became non-finite at iteration 0 (simulated)")
        return real(Y, cfg)

    monkeypatch.setattr(fit_mod, "_fit_once", flaky)
    with pytest.warns(RuntimeWarning, match="seed-family"):
        res = fit_cebmf(Y, FitConfig(K=4, maxit=30, internal_epoch=10, seed=1))
    fam = res.guard_diagnostics["seed_family"]
    assert fam[0]["diverged"] is True
    assert not fam[1]["diverged"]
    assert res.guard_diagnostics["seed_applied"] == 1 + fit_mod.SEED_FAMILY_STRIDE


def test_clean_fit_no_fallback_no_warning():
    """A fit that never diverges: no RuntimeWarning, single-attempt log."""
    counts, _, _, _ = make_counts(n_cells=80, n_genes=60, n_programs=2, seed=1)
    Y = torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        res = fit_cebmf(Y, FitConfig(K=4, maxit=30, internal_epoch=10, seed=1))
    fam = res.guard_diagnostics["seed_family"]
    assert len(fam) == 1 and fam[0]["diverged"] is False
    assert res.guard_diagnostics["seed_applied"] == 1
    assert np.isfinite(res.objective[-1])


def test_fixed_s_known_noise_runs():
    """FitConfig.S fixes the residual SD (KNOWN noise) and the fit stays finite."""
    counts, _, _, _ = make_counts(n_cells=80, n_genes=60, n_programs=2, seed=3)
    Y = torch.tensor(ShiftedLog().fit_transform(counts), dtype=torch.float32)
    res = fit_cebmf(Y, FitConfig(K=4, maxit=30, internal_epoch=10, seed=3, S=0.1))
    assert np.isfinite(res.objective[-1])
    assert res.effective_k >= 1
    assert np.all(np.isfinite(res.tau)) and np.all(res.tau <= fit_mod.TAU_CAP)
