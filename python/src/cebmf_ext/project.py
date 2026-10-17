"""Chunked single-cell projection with fixed gene loadings (M2).

Given a discovery fit's F (genes x K) and tau, infer cell scores L for a
sparse count matrix in row-chunks, never materializing the full dense matrix.

Exactness map (design doc 03 §3, as amended by M2 validation):
- The likelihood factorizes over rows given (F, tau) -> chunked inference is
  exact for the conditional L-update.
- The package's L-update REFITS the L-prior on the current column of
  statistics each call (cebmf.py:347-360). Under chunking the prior EM sees
  the chunk's values instead of all N -> the chunked path is approximate with
  respect to the full-batch prior fit. Agreement is measured empirically
  (test_chunk_equivalence), not assumed.
- tau is taken from the discovery fit and held fixed during projection.

Implementation reuses the package's own coordinate updates: run iter_once on
a chunk-scoped model, then restore the frozen F/F2/tau (the F-update and
tau-update results are discarded each sweep).
"""
from __future__ import annotations

import numpy as np
import torch

from cebmf_torch import cEBMF
from cebmf_torch.cebmf.cebmf import NoiseType

from .io import dense_chunk


def _conditional_L_sweeps(Y_chunk: torch.Tensor, F: torch.Tensor, tau: float,
                          K: int, n_sweeps: int, internal_epoch: int = 1) -> torch.Tensor:
    """L-updates with F and tau frozen, using the package's own update path."""
    model = cEBMF(
        Y_chunk, K=K, prior_L="laplace", prior_F="laplace",
        noise_type=NoiseType.CONSTANT,  # tau frozen below anyway; CONSTANT avoids per-chunk tau drift
        allow_backfitting=False,
    )
    # deterministic zero init for L, then freeze F at the discovery solution
    model.initialise_factors(
        L=torch.zeros(Y_chunk.shape[0], K, dtype=torch.float32),
        F=F.clone(),
    )
    for _ in range(n_sweeps):
        model.iter_once()
        # discard the F-update and tau-update; keep only the L-update effect
        model.F = F.clone()
        model.F2 = (F * F).clone()
        model.tau = torch.tensor(float(tau))
        model.tau_map = torch.full(Y_chunk.shape, float(tau))
        model._recompute_residual()
    return model.L.clone()


def project_scores(
    X: "scipy.sparse.csr_matrix",
    F: np.ndarray,
    tau: float,
    n_sweeps: int = 5,
    chunk_size: int = 10000,
    transform=None,
) -> np.ndarray:
    """Project gene loadings onto cells in chunks.

    X: cells x genes SPARSE counts. F: genes x K (discovery loadings, in
    transformed space). tau: scalar precision from the discovery fit
    (COLUMN_WISE tau is per-gene; the frozen-tau path here uses the scalar
    mean precision — per-gene tau_map support lands with the M2 benchmark).
    transform: optional cebmf_ext.preproc.ShiftedLog with fitted state; each
    chunk is transformed with its own library-size size factors (the recipe
    is per-cell and needs no cross-cell state).
    Returns L (cells x K) dense.
    """
    F_t = torch.tensor(F, dtype=torch.float32)
    K = F.shape[1]
    n = X.shape[0]
    out = np.zeros((n, K), dtype=np.float32)
    for start in range(0, n, chunk_size):
        sl = slice(start, min(start + chunk_size, n))
        Xc = dense_chunk(X, sl)
        if transform is not None:
            Xc = transform.transform(Xc)
        Yc = torch.tensor(np.asarray(Xc), dtype=torch.float32)
        out[sl] = _conditional_L_sweeps(Yc, F_t, tau, K, n_sweeps).numpy()
    return out
