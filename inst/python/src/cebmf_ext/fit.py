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
from cebmf_torch.priors.point import PointBuilder

# --- R1 numerical guard set (M6 revision R1; ledger deviations 6-7) ---------
# Registered values, re-integrated 2026-10-18 from the audit ledger:
# (a) sebetahat floor 1e-4 before every point-prior (laplace/exp/gbinary) fit.
#     Mechanism (GZMKpos root cause, ledger 2026-10-18): single-donor-spike
#     genes are exactly fittable -> residual variance -> 0 -> sebetahat -> 0
#     -> NaN inside the laplace EBNM optimizer. Class-level guard on
#     PointBuilder.fit, the same wrapper layer as the tau cap.
# (b) tau precision cap 1e8 on the returned precision (the backend's own cap
#     is 1e12; the wrapper clamps tighter so downstream precision-weighted
#     code cannot overflow).
# (c) explicit seed application: np.random.seed + torch.manual_seed.
#     (FitConfig.seed was previously a no-op.)
# (d) seed-family fallbacks: on objective divergence, retry with a seed family
#     and record diagnostics (registered deviation).
SEBETAHAT_FLOOR = 1e-4
TAU_CAP = 1e8
SEED_FAMILY_ATTEMPTS = 3
SEED_FAMILY_STRIDE = 1000


def _apply_numerical_guards() -> None:
    """Install the class-level sebetahat floor on PointBuilder.fit (idempotent)."""
    if getattr(PointBuilder, "_cellprograms_sebetahat_guard", False):
        return
    original_fit = PointBuilder.fit

    def guarded_fit(self, X=None, betahat=None, sebetahat=None, model_param=None,
                    internal_epoch=None, device=None):
        if sebetahat is not None:
            sebetahat = torch.clamp(sebetahat, min=SEBETAHAT_FLOOR)
        return original_fit(self, X=X, betahat=betahat, sebetahat=sebetahat,
                            model_param=model_param,
                            internal_epoch=internal_epoch, device=device)

    PointBuilder.fit = guarded_fit
    PointBuilder._cellprograms_sebetahat_guard = True


_apply_numerical_guards()

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
    # Fixed residual SD (var_floor exposure, ledger 2026-10-18). When set,
    # cebmf_torch forces NoiseType.KNOWN: the noise variance is FIXED at S**2
    # and `noise` is ignored (cebmf.py constructor contract). Scalar only from
    # the R bridge; tensor S stays a python-side option.
    S: float | None = None

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
    guard_diagnostics: dict = field(default_factory=dict)  # R1 deviation record


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


def _fit_once(Y: torch.Tensor, cfg: FitConfig) -> FitResult:
    """Single guarded fit attempt. Raises RuntimeError on objective divergence."""
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
        S=cfg.S,
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
    tau = fit.tau.cpu().numpy() if torch.is_tensor(fit.tau) else np.asarray(fit.tau)
    # R1(b): tau precision cap — see module docstring guard notes.
    tau_capped = np.minimum(tau, TAU_CAP)
    return FitResult(
        L=L,
        F=fit.F.cpu().numpy(),
        tau=tau_capped,
        objective=trace,
        direction=_objective_direction(trace),
        effective_k=effective_k(L),
        converged=converged,
        runtime_s=time.time() - t0,
        peak_mem_gb=peak_gb,
        config=cfg,
        guard_diagnostics={
            "tau_cap_applied": bool(np.any(tau > TAU_CAP)),
            "tau_max_raw": float(np.max(tau)) if tau.size else None,
        },
    )


def fit_cebmf(Y: torch.Tensor, cfg: FitConfig | None = None) -> FitResult:
    """Fit cEBMF with ELBO-plateau stopping. Y: obs x genes float tensor, NaN = missing.

    R1 guard set applied here (registered deviations, ledger 2026-10-18):
    explicit seed application, sebetahat floor (class-level, module import),
    tau cap 1e8, and a seed-family fallback on objective divergence. Divergence
    after the seed family is exhausted re-raises with the attempt diagnostics.
    """
    cfg = cfg or FitConfig()
    attempts: list[dict] = []
    result: FitResult | None = None
    last_err: RuntimeError | None = None
    for attempt in range(SEED_FAMILY_ATTEMPTS):
        seed_k = cfg.seed + attempt * SEED_FAMILY_STRIDE
        # R1(c): explicit seed application. FitConfig.seed was a no-op before
        # this guard; both RNG streams are seeded for reproducibility.
        np.random.seed(seed_k)
        torch.manual_seed(seed_k)
        try:
            result = _fit_once(Y, cfg)
        except RuntimeError as err:
            attempts.append({"attempt": attempt, "seed": seed_k,
                             "diverged": True, "error": str(err)[:200]})
            last_err = err
            continue
        attempts.append({"attempt": attempt, "seed": seed_k, "diverged": False})
        result.guard_diagnostics["seed_family"] = attempts
        result.guard_diagnostics["seed_applied"] = seed_k
        if attempt > 0:
            # Registered deviation: a fallback occurred; surface it loudly.
            import warnings

            warnings.warn(
                f"cEBMF diverged on seed {cfg.seed}; recovered with seed-family "
                f"fallback (seed {seed_k}, attempt {attempt}).",
                RuntimeWarning,
                stacklevel=2,
            )
        return result
    raise RuntimeError(
        f"cEBMF diverged for all {SEED_FAMILY_ATTEMPTS} seed-family attempts "
        f"(seeds {cfg.seed}..{seed_k}); last error: {last_err}"
    ) from last_err
