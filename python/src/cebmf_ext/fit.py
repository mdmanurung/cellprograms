"""Thin cEBMF fitting wrapper: ELBO-plateau stopping, effective K, resource capture.

The wrapper never changes the cebmf_torch objective; it adds the convergence
control the package lacks (fixed-maxit only, audit 01) at the adapter level by
calling fit() in internal-epoch chunks and watching the recorded objective.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import torch

from cebmf_torch import cEBMF
from cebmf_torch.cebmf.cebmf import NoiseType

NOISE_MAP = {
    "constant": NoiseType.CONSTANT,
    "row_wise": NoiseType.ROW_WISE,
    "column_wise": NoiseType.COLUMN_WISE,
}


@dataclass
class FitConfig:
    K: int = 10
    prior_L: str = "laplace"
    prior_F: str = "laplace"
    noise: str = "column_wise"
    allow_backfitting: bool = True
    prune_thresh: float = 1.0 - 1e-3
    internal_epoch: int = 10
    maxit: int = 100  # total iteration budget across chunks
    elbo_tol: float = 1e-5  # relative objective improvement considered a plateau
    elbo_patience: int = 3  # consecutive plateau checks before stopping
    seed: int = 0
    init_L: np.ndarray | None = None  # obs x K user init (M3 refinement hook)
    init_F: np.ndarray | None = None  # genes x K user init

    def to_dict(self) -> dict:
        return dict(self.__dict__)


@dataclass
class FitResult:
    L: np.ndarray  # obs x K
    F: np.ndarray  # genes x K
    tau: np.ndarray
    objective: list  # objective trace (LOSS = negative ELBO; source-verified)
    direction: str  # always "loss" for cebmf_torch (cebmf.py:427-428)
    effective_k: int
    converged: bool  # True if stopped by plateau, False if maxit exhausted
    runtime_s: float
    peak_mem_gb: float
    config: FitConfig = field(repr=False, default=None)


def _objective_direction(trace: list[float]) -> str:
    """cebmf_torch records LOSS = -loglik + KL and minimizes it
    (cebmf.py:427-428, source-verified). Endpoint-based auto-detection is
    unreliable because the KL terms start at zero (priors unfitted), making
    the first recorded values artificially low — see tests/test_elbo_monotone."""
    return "loss"


def effective_k(L: np.ndarray, tol: float = 1e-8) -> int:
    """Surviving factors: columns with non-degenerate norm after pruning."""
    norms = np.linalg.norm(L, axis=0)
    return int((norms > tol * max(norms.max(), 1e-300)).sum())


def fit_cebmf(Y: torch.Tensor, cfg: FitConfig | None = None) -> FitResult:
    """Fit cEBMF with ELBO-plateau stopping. Y: obs x genes float tensor, NaN = missing."""
    cfg = cfg or FitConfig()
    if Y.dtype != torch.float32:
        Y = Y.float()
    t0 = time.time()
    model = cEBMF(
        Y,
        K=cfg.K,
        prior_L=cfg.prior_L,
        prior_F=cfg.prior_F,
        internal_epoch=cfg.internal_epoch,
        allow_backfitting=cfg.allow_backfitting,
        prune_thresh=cfg.prune_thresh,
        noise_type=NOISE_MAP[cfg.noise],
    )
    if cfg.init_L is not None and cfg.init_F is not None:
        model.initialise_factors(
            L=torch.as_tensor(cfg.init_L, dtype=torch.float32, device=model.device),
            F=torch.as_tensor(cfg.init_F, dtype=torch.float32, device=model.device),
        )
    trace: list[float] = []
    converged = False
    it = 0
    while it < cfg.maxit:
        step = min(cfg.internal_epoch, cfg.maxit - it)
        fit = model.fit(maxit=step)
        new = [float(x) for x in fit.history_obj] if fit.history_obj is not None else []
        # history_obj is cumulative (fit() does not reset obj; only init and
        # pruning events do — cebmf.py:520,650), so keep the longest trace.
        trace = new if len(new) >= len(trace) else trace + new
        # Fail loudly on objective divergence: cebmf_torch can diverge silently
        # (loss -> NaN) on some inputs (observed: small synthetic, over-fitted K;
        # ledger 2026-10-17). Silent NaN factors are exactly the failure mode
        # this audit program exists to prevent.
        if trace and not np.isfinite(trace[-1]):
            raise RuntimeError(
                f"cEBMF objective became non-finite at iteration {it} "
                f"(config: K={cfg.K}, priors={cfg.prior_L}/{cfg.prior_F}, "
                f"noise={cfg.noise}). This is a backend divergence, not a "
                "missing-data effect; reduce K, check input scaling, or report upstream."
            )
        it += step
        if len(trace) >= 2:
            scale = max(abs(trace[-1]), abs(trace[0]), 1e-300)
            gain = (trace[-2] - trace[-1]) / scale  # loss decrease = improvement
            if gain < cfg.elbo_tol:
                tail = trace[-(cfg.elbo_patience + 1):]
                if len(tail) == cfg.elbo_patience + 1:
                    deltas = (np.diff(tail) * -1) / scale  # improvements
                    if np.all(deltas < cfg.elbo_tol):
                        converged = True
                        break
    import resource

    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # Linux: KB
    L = fit.L.cpu().numpy()
    return FitResult(
        L=L,
        F=fit.F.cpu().numpy(),
        tau=fit.tau.cpu().numpy() if torch.is_tensor(fit.tau) else np.asarray(fit.tau),
        objective=trace,
        direction=_objective_direction(trace),
        effective_k=effective_k(L),
        converged=converged,
        runtime_s=time.time() - t0,
        peak_mem_gb=peak_gb,
        config=cfg,
    )
