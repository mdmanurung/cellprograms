"""M3 synthetic generator (pre-registered: docs/cebmf-extension/04-experiment-plan.md §1).

5,000 cells x 1,500 genes; 50 donors x 100 cells. Truth: 3 broad programs (150 markers
each, block W; donor-level usage with donor random effect sd 0.3) + 1 rare program
(5% of cells, 80 markers) + 2 null usage columns with NO loading support (unobservable;
they test over-fitted-K pruning). Counts ~ Poisson(gene-specific abundance x exp(program
effects)) — UMI-style heteroscedastic (Amendment A3); nested Bernoulli dropout
(u drawn once per rep, dropout iff u < rate). Donor covariate z_age; optional planted
moderation on program-1 usage amplitude (beta; E-M3.1 runs beta=0).

Deterministic per rep: every (rep, rate, beta) regenerates identically from seed
100000+rep — no intermediate data files.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

N_CELLS, N_GENES, N_DONORS, CPD = 5000, 1500, 50, 100
K_BROAD, K_NULL = 3, 2
BROAD_BLOCK, RARE_BLOCK = 150, 80
N_RARE = N_CELLS // 20  # 250 cells = 5%
N_TEST_DONORS = 10


def generate(rep: int, rate: float, beta: float = 0.0) -> dict:
    rng = np.random.default_rng(100000 + rep)
    donor_ids = np.repeat(np.arange(N_DONORS), CPD)
    z_age = rng.standard_normal(N_DONORS)
    donor_eff = rng.standard_normal((N_DONORS, K_BROAD)) * 0.3
    z = rng.standard_normal((N_CELLS, K_BROAD)) + donor_eff[donor_ids]
    if beta:
        z[:, 0] = z[:, 0] * (1.0 + beta * z_age[donor_ids])
    W = np.zeros((N_GENES, K_BROAD))
    for k in range(K_BROAD):
        W[k * BROAD_BLOCK:(k + 1) * BROAD_BLOCK, k] = 1.0
    broad_markers = [list(range(k * BROAD_BLOCK, (k + 1) * BROAD_BLOCK)) for k in range(K_BROAD)]
    rare_markers = list(range(K_BROAD * BROAD_BLOCK, K_BROAD * BROAD_BLOCK + RARE_BLOCK))
    rare_cells = np.sort(rng.choice(N_CELLS, size=N_RARE, replace=False))
    z_rare = np.zeros(N_CELLS)
    z_rare[rare_cells] = rng.standard_normal(N_RARE) * 1.5 + 2.0
    z_null = rng.standard_normal((N_CELLS, K_NULL))  # no loading support

    # Amendment A3 (pre-run): Poisson UMI-style counts with gene-specific abundance,
    # replacing the E-style exp(Gaussian) lognormal model. The pilot showed the
    # homoscedastic lognormal generator violates the registered CP10K+log1p+scale
    # preprocessing (gene-wise scaling amplifies 1,450 uniform noise genes to
    # variance 1, drowning the marker blocks; SVD capped at 0.66-0.80 regardless of
    # SNR). Real UMI counts are heteroscedastic (log-scale variance ~ 1/mean), which
    # is the structure the recipe expects. Usage noise 0.2 (A2) retained; the base
    # Gaussian term is REPLACED by Poisson sampling. Gates/arms/metrics unchanged.
    log_mu = rng.normal(3.0, 1.0, N_GENES)  # gene abundance, median ~20 counts/cell
    lam = np.exp(log_mu)[None, :] * np.exp(
        (z + rng.standard_normal((N_CELLS, K_BROAD)) * 0.2) @ W.T
        + (z_rare + rng.standard_normal(N_CELLS) * 0.2)[:, None] * (np.isin(np.arange(N_GENES), rare_markers))
    )
    # Decouple library size from program usage (A3): without this, usage drives
    # total counts and CP10K turns it into a dominant compositional-dilution
    # direction (pilot: PC1 = lib size, cor 0.84, sv 2x next; SVD ceiling 0.69-0.83).
    # Real lib-size variation reflects cell size/RNA content, not planted programs.
    row_tot = lam.sum(axis=1, keepdims=True)
    lam = lam * (row_tot.mean() / row_tot)
    counts = rng.poisson(lam).astype(np.float32)

    u = rng.random((N_CELLS, N_GENES))  # nested dropout masks (E design)
    counts[u < rate] = 0.0

    test_donors = np.sort(rng.choice(N_DONORS, N_TEST_DONORS, replace=False))
    test_mask = np.isin(donor_ids, test_donors)

    Z_true = np.column_stack([z, z_rare])  # 4 truth columns: P1,P2,P3 broad + P4 rare
    markers = broad_markers + [rare_markers]
    return dict(
        counts=sparse.csr_matrix(counts, dtype=np.float32),
        Z_true=Z_true, markers=markers, donor_ids=donor_ids, z_age=z_age,
        test_mask=test_mask, rare_cells=rare_cells, rate=rate, beta=beta, rep=rep,
    )
