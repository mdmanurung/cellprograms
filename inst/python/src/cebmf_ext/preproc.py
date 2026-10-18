"""Replayable shifted-log preprocessing (flashier single-cell vignette recipe).

y = log1p(x / (alpha * s_i)) with library-size size factors s_i (mean 1),
following Ahlmann-Eltz & Huber 2023 as used by the flashier single-cell
vignette (pinned fd44811). Gene centering/scaling is computed once on the
discovery matrix and replayed exactly on projection chunks.

The recipe object is JSON-serializable so the result contract can replay it.
HVG selection is external: the adapter consumes a pre-selected gene index
list (the granularity builders own HVG choice).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np


@dataclass
class ShiftedLog:
    """Shifted-log transform with stored gene center/scale for replay."""

    alpha: float = 1.0
    # Fitted state (filled by fit_transform on the discovery matrix).
    gene_mean: np.ndarray | None = None
    gene_sd: np.ndarray | None = None
    size_factors: np.ndarray | None = None  # discovery-matrix size factors
    n_hvg: int = -1

    def _log_matrix(self, counts: np.ndarray, sf: np.ndarray) -> np.ndarray:
        if (counts < 0).any():
            raise ValueError("counts must be non-negative")
        return np.log1p(counts / (self.alpha * sf[:, None]))

    def fit_transform(self, counts: np.ndarray) -> np.ndarray:
        """Fit on the discovery matrix (cells x genes, counts) and return Y."""
        lib = counts.sum(axis=1).astype(float)
        if (lib <= 0).any():
            raise ValueError(
                f"{int((lib <= 0).sum())} cells have zero counts over the gene panel; "
                "drop them before fitting (structural zeros cannot be normalized)"
            )
        sf = lib / lib.mean()
        Y = self._log_matrix(counts, sf)
        self.gene_mean = Y.mean(axis=0)
        self.gene_sd = Y.std(axis=0, ddof=0)
        self.gene_sd[self.gene_sd == 0] = 1.0  # constant genes: center only
        self.size_factors = sf
        self.n_hvg = counts.shape[1]
        return self._scale(Y)

    def transform(self, counts: np.ndarray, size_factors: np.ndarray | None = None) -> np.ndarray:
        """Replay on a chunk. size_factors: per-cell factors for the chunk
        (recomputed from the chunk's own library size by default — the
        vignette recipe is per-cell and needs no cross-cell state)."""
        if self.gene_mean is None:
            raise RuntimeError("call fit_transform before transform")
        if size_factors is None:
            lib = counts.sum(axis=1).astype(float)
            if (lib <= 0).any():
                raise ValueError("chunk contains zero-library-size cells")
            size_factors = lib / lib.mean()
        return self._scale(self._log_matrix(counts, size_factors))

    def _scale(self, Y: np.ndarray) -> np.ndarray:
        return (Y - self.gene_mean) / self.gene_sd

    def to_dict(self) -> dict:
        return {
            "recipe": "shifted_log",
            "alpha": self.alpha,
            "size_factor": "library_size_mean1",
            "gene_center": True,
            "gene_scale": True,
            "n_hvg": int(self.n_hvg),
            "gene_mean": None if self.gene_mean is None else self.gene_mean.tolist(),
            "gene_sd": None if self.gene_sd is None else self.gene_sd.tolist(),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ShiftedLog":
        obj = cls(alpha=d["alpha"])
        obj.n_hvg = d["n_hvg"]
        obj.gene_mean = None if d["gene_mean"] is None else np.asarray(d["gene_mean"])
        obj.gene_sd = None if d["gene_sd"] is None else np.asarray(d["gene_sd"])
        return obj

    def dumps(self) -> str:
        return json.dumps(self.to_dict())
