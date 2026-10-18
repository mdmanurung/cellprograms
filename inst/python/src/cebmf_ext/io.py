"""Sparse loaders: never materialize the full dense matrix at load time."""
from __future__ import annotations

import numpy as np
import scipy.io as sio
from scipy.sparse import csr_matrix


def load_mtx(path: str) -> csr_matrix:
    """Load a MatrixMarket counts matrix as CSR (cells x genes)."""
    X = sio.mmread(path)
    if not hasattr(X, "tocsr"):
        X = csr_matrix(X)
    return X.tocsr()


def load_lines(path: str) -> list[str]:
    with open(path) as f:
        return [line.rstrip("\n") for line in f]


def dense_chunk(X: csr_matrix, rows: slice | np.ndarray) -> np.ndarray:
    """Materialize ONE chunk densely (the only dense materialization allowed)."""
    return np.asarray(X[rows].todense())
