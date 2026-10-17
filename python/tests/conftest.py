"""Shared synthetic-counts generator for the cebmf_ext test suite."""
from __future__ import annotations

import numpy as np


def make_counts(n_cells=150, n_genes=40, n_programs=3, seed=0, dropout_rate=0.0):
    """Block-structured counts with planted programs.

    Returns (counts, W_true genes x p, S_true cells x p, program_genes list).
    """
    rng = np.random.default_rng(seed)
    genes_per_program = max(3, n_genes // (2 * n_programs))
    W = np.zeros((n_genes, n_programs))
    program_genes = []
    for p in range(n_programs):
        idx = rng.choice(n_genes, size=genes_per_program, replace=False)
        W[idx, p] = rng.uniform(0.5, 1.5, size=genes_per_program)
        program_genes.append(idx)
    S = np.zeros((n_cells, n_programs))
    for p in range(n_programs):
        S[rng.choice(n_cells, size=n_cells // 3, replace=False), p] = rng.uniform(1.0, 2.5, size=n_cells // 3)
    lam = np.exp(np.clip(S @ W.T, 0.0, 4.0)) + 0.05
    counts = rng.poisson(lam).astype(float)
    if dropout_rate > 0:
        mask = rng.random(counts.shape) < dropout_rate
        counts[mask] = 0.0
    return counts, W, S, program_genes
