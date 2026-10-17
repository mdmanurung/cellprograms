"""Normalization variants for Amendment A4 (pre-registered, 04d-amendment-A4).

Alternatives to the registered shifted-log recipe (control) for the S1
(single-cell) and S3 (metacell) strategy arms:

- PFLog  (Booeshaghi et al. 2022, "shifted CLR", Anscombe-calibrated):
  y = rowcenter( log1p(4 * alpha_hat * x / s_c) ),  s_c = library size / mean.
  alpha_hat estimated by pooled NB moment matching on depth-normalized TRAIN
  counts (var = mu + alpha*mu^2), median over expressed genes, floored at 0.25
  (Poisson-limit fallback keeps 4*alpha_hat >= 1, i.e. the conventional
  pseudocount-1 scale). alpha_hat is recorded per run.
- CLRLog (plain CLR on the log1p scale): y = rowcenter( log1p(x) ) — the
  per-cell geometric-mean centering IS the depth normalization; no size
  factors.

All three arms share the registered gene center+scale stage (fit on TRAIN,
replayed on test chunks), so arms differ only in the count-scale/centering
geometry. Interfaces match cebmf_ext.preproc.ShiftedLog (fit_transform /
transform), so project_scores consumes them unchanged.
"""
from __future__ import annotations

import numpy as np

ALPHA_FLOOR = 0.25


def _row_center(Y: np.ndarray) -> np.ndarray:
    return Y - Y.mean(axis=1, keepdims=True)


def _lib_size_factors(counts: np.ndarray) -> np.ndarray:
    lib = counts.sum(axis=1).astype(float)
    if (lib <= 0).any():
        raise ValueError(
            f"{int((lib <= 0).sum())} rows have zero counts over the gene panel; "
            "drop them before fitting (structural zeros cannot be normalized)"
        )
    return lib / lib.mean()


def estimate_alpha(counts: np.ndarray) -> float:
    """Pooled NB overdispersion by moment matching on depth-normalized counts."""
    sf = _lib_size_factors(counts)
    X = counts / sf[:, None]
    mu = X.mean(axis=0)
    var = X.var(axis=0, ddof=0)
    expressed = mu > 0.01 * mu.max()
    if not expressed.any():
        return ALPHA_FLOOR
    mu_e, var_e = mu[expressed], var[expressed]
    with np.errstate(divide="ignore", invalid="ignore"):
        a_j = (var_e - mu_e) / (mu_e ** 2)
    a_j = a_j[np.isfinite(a_j) & (a_j > 0)]
    if a_j.size == 0:
        return ALPHA_FLOOR
    return float(max(np.median(a_j), ALPHA_FLOOR))


class _GeneScaleMixin:
    """Registered gene center+scale stage (identical to ShiftedLog)."""

    def _fit_gene_stage(self, Y: np.ndarray) -> np.ndarray:
        self.gene_mean = Y.mean(axis=0)
        self.gene_sd = Y.std(axis=0, ddof=0)
        self.gene_sd[self.gene_sd == 0] = 1.0
        return (Y - self.gene_mean) / self.gene_sd

    def _apply_gene_stage(self, Y: np.ndarray) -> np.ndarray:
        if self.gene_mean is None:
            raise RuntimeError("call fit_transform before transform")
        return (Y - self.gene_mean) / self.gene_sd


class PFLog(_GeneScaleMixin):
    """Shifted CLR with Anscombe-calibrated pseudocount (Booeshaghi 2022)."""

    def __init__(self) -> None:
        self.alpha_hat: float | None = None
        self.gene_mean: np.ndarray | None = None
        self.gene_sd: np.ndarray | None = None

    def _core(self, counts: np.ndarray) -> np.ndarray:
        sf = _lib_size_factors(counts)
        return _row_center(np.log1p(4.0 * self.alpha_hat * counts / sf[:, None]))

    def fit_transform(self, counts: np.ndarray) -> np.ndarray:
        self.alpha_hat = estimate_alpha(counts)
        return self._fit_gene_stage(self._core(counts))

    def transform(self, counts: np.ndarray, size_factors: np.ndarray | None = None) -> np.ndarray:
        return self._apply_gene_stage(self._core(counts))

    def to_dict(self) -> dict:
        return {"recipe": "pflog_shifted_clr", "alpha_hat": self.alpha_hat,
                "alpha_floor": ALPHA_FLOOR, "row_center": True,
                "gene_center": True, "gene_scale": True}


class CLRLog(_GeneScaleMixin):
    """Plain CLR on the log1p scale (per-cell geometric-mean centering)."""

    def __init__(self) -> None:
        self.gene_mean: np.ndarray | None = None
        self.gene_sd: np.ndarray | None = None

    @staticmethod
    def _core(counts: np.ndarray) -> np.ndarray:
        if (counts < 0).any():
            raise ValueError("counts must be non-negative")
        return _row_center(np.log1p(counts))

    def fit_transform(self, counts: np.ndarray) -> np.ndarray:
        return self._fit_gene_stage(self._core(counts))

    def transform(self, counts: np.ndarray, size_factors: np.ndarray | None = None) -> np.ndarray:
        return self._apply_gene_stage(self._core(counts))

    def to_dict(self) -> dict:
        return {"recipe": "clr_log1p", "row_center": True,
                "gene_center": True, "gene_scale": True}


def make_recipe(norm: str):
    """Factory used by run_arm.py --norm."""
    if norm == "shiftedlog":
        from cebmf_ext.preproc import ShiftedLog
        return ShiftedLog()
    if norm == "pflog":
        return PFLog()
    if norm == "clr":
        return CLRLog()
    raise ValueError(f"unknown norm: {norm}")
