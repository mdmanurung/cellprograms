"""Test 5 (plan, as amended): chunked vs dense conditional-L agreement.

AMENDMENT (recorded before the M2 benchmark): the plan's "< 1e-4 relative"
criterion assumed the prior state is frozen during projection. The package's
L-update refits the L-prior on the current statistics each call
(cebmf.py:347-360), so under chunking the prior EM sees chunk-level values —
the chunked path is approximate with respect to the full-batch prior fit.
The criterion is therefore: chunked and dense conditional-L sweeps converge
to the same conditional optimum with measured agreement reported; the test
criterion is EMPIRICAL with the mechanism characterized (2026-10-17,
synthetic 1000x60, 3 planted programs): signal-factor shape error is
prior-refit-driven and CONVERGED (stable across 5-40 sweeps), scaling with
chunk count — measured 3.0-10.9% at 10 chunks, 1.9-6.4% at 4 chunks,
1.3-4.7% at 2 chunks, 0% at 1 chunk (exact sanity check). Test criterion:
matched |cor| >= 0.99 for >= 3 signal factors and relative error <= 0.15
at 4 chunks. Upstream feature request: posterior-only updates with frozen
prior state (would restore near-exactness).
"""
import numpy as np
import torch
from scipy import sparse

from cebmf_ext import FitConfig, fit_cebmf
from cebmf_ext.preproc import ShiftedLog
from cebmf_ext.project import _conditional_L_sweeps, project_scores

from conftest import make_counts


def test_chunk_equivalence():
    counts, _, S_true, _ = make_counts(n_cells=1000, n_genes=60, n_programs=3, seed=7)
    recipe = ShiftedLog()
    Y = recipe.fit_transform(counts)

    # discovery fit (dense) to get F and tau
    res = fit_cebmf(torch.tensor(Y, dtype=torch.float32),
                    FitConfig(K=5, maxit=40, internal_epoch=10, seed=7))
    F = res.F
    tau = float(np.mean(1.0 / np.asarray(res.tau))) if np.ndim(res.tau) else float(res.tau)

    # dense conditional-L sweep (single chunk)
    L_dense = _conditional_L_sweeps(
        torch.tensor(Y, dtype=torch.float32), torch.tensor(F), tau, K=F.shape[1], n_sweeps=5
    ).numpy()

    # chunked projection (4 chunks) from the SPARSE counts, replaying preprocessing
    X_sp = sparse.csr_matrix(counts)
    L_chunked = project_scores(X_sp, F, tau, n_sweeps=5, chunk_size=250, transform=recipe)

    # matched |cor| per factor — the criterion applies to SIGNAL factors
    # (planted programs); excess noise factors (K=5 > 3 planted) are
    # prior-refit-dominated and not required to be chunk-stable (same
    # phenomenon as the masking test).
    A = L_dense - L_dense.mean(0)
    B = L_chunked - L_chunked.mean(0)
    C = np.corrcoef(A.T, B.T)[: F.shape[1], F.shape[1]:]
    C = np.abs(C)
    matched = C.max(axis=1)
    n_signal = int((matched >= 0.99).sum())
    assert n_signal >= 3, f"only {n_signal}/3 signal factors chunk-stable: {matched}"

    # relative error after centering + sign/scale alignment (signal factors
    # only). Centering matches the correlation criterion's invariance; the
    # per-chunk prior refit may shift column means without changing shape.
    for k in np.where(matched >= 0.99)[0]:
        j = int(np.argmax(C[k]))
        a = L_dense[:, k] - L_dense[:, k].mean()
        b = L_chunked[:, j] - L_chunked[:, j].mean()
        scale = np.dot(a, b) / max(np.dot(b, b), 1e-300)
        rel = np.linalg.norm(a - scale * b) / max(np.linalg.norm(a), 1e-300)
        assert rel <= 0.15, f"factor {k}: relative error {rel:.3e} > 0.15"
