"""M3 metrics: Hungarian matching, top-50 Jaccard, held-out prediction."""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def match_cor_matrix(L_est: np.ndarray, Z_true: np.ndarray) -> np.ndarray:
    """|Pearson cor| between every estimated score column and truth column."""
    A = L_est - L_est.mean(0)
    B = Z_true - Z_true.mean(0)
    C = np.corrcoef(A.T, B.T)[: A.shape[1], A.shape[1]:]
    return np.abs(np.nan_to_num(C, nan=0.0))


def hungarian_match(L_est: np.ndarray, Z_true: np.ndarray) -> list[dict]:
    """Optimal 1-1 match on |cor|; returns one row per truth column.

    C is (K_est x K_true); linear_sum_assignment returns (est_rows, truth_cols).
    """
    C = match_cor_matrix(L_est, Z_true)
    est_rows, truth_cols = linear_sum_assignment(-C)
    out = []
    for e, t in zip(est_rows, truth_cols):
        out.append(dict(truth_col=int(t), est_col=int(e), match_cor=float(C[e, t])))
    out.sort(key=lambda d: d["truth_col"])
    return out


def top_jaccard(F_est: np.ndarray, est_col: int, marker_set: set, top: int = 50) -> float:
    top_genes = set(np.argsort(-np.abs(F_est[:, est_col]))[:top])
    inter = len(top_genes & marker_set)
    return inter / max(len(top_genes | marker_set), 1)


def heldout_prediction(L_test: np.ndarray, F: np.ndarray, Y_test: np.ndarray,
                       mask: np.ndarray) -> dict:
    pred = L_test @ F.T
    err = pred[mask] - Y_test[mask]
    rmse = float(np.sqrt(np.mean(err ** 2)))
    p = float(np.corrcoef(pred[mask], Y_test[mask])[0, 1]) if mask.sum() > 2 else float("nan")
    return dict(heldout_pearson=p, heldout_rmse=rmse, n_masked=int(mask.sum()))


def fixed_entry_mask(shape: tuple[int, int], rep: int, frac: float = 0.1) -> np.ndarray:
    """Fixed random 10% entry mask for held-out prediction (seeded per rep)."""
    rng = np.random.default_rng(500000 + rep)
    m = rng.random(shape) < frac
    m[0, 0] = True  # guard against an all-False draw
    return m
