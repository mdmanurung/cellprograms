# 04c — Amendment A3 (2026-10-17, pre-run): Poisson UMI generator + ceiling-relative gates

Registered BEFORE any M3 battery result exists (S1/S2 pilots only). Supersedes the
uncommitted A2 generator change; A2's committed diagnosis is retained as history.

## Pilot findings (rep 1, rate 0, beta 0)

1. The Hungarian matching in m3_metrics.py had its assignment indices swapped
   (est/truth labels transposed; Jaccard looked up the wrong F column). Fixed; pilot
   metric LABELS were invalid, correlation values valid.
2. Pseudobulk aggregation averaged over all 50 donor ids including the 10 held-out
   donors (10 all-NaN rows). Fixed: aggregate over training donors only.
3. Generator diagnosis (deeper than A2): with the E-style lognormal generator, SVD on
   the registered recipe caps at 0.66-0.80 per planted program regardless of SNR
   (A2's sd fix changed nothing). Measured mechanism, in order: (a) PC1 was a LIBRARY-
   SIZE direction (cor 0.84, singular value 2x next) because program usage drove total
   counts — fixed by decoupling library size from usage (row-normalized lambda);
   (b) the remaining ceiling is COMPOSITIONAL COUPLING between co-occurring programs:
   the relative-abundance operator couples the blocks' observable signals through the
   shared log-denominator, so top SVD PCs are mixtures (e.g. PC2 = 0.69 P2 + 0.70 P3).
   Per-gene signal is healthy (0.87-0.90 broad, 0.71 rare) and within-block PCs reach
   0.94 — the loss is inherent to multi-program cells under compositional
   normalization, which is REALISTIC (real multi-program cells behave the same way).

## Changes

1. **Generator**: Poisson UMI-style counts — counts ~ Poisson(lambda), lambda =
   gene-specific abundance (log mu ~ N(3, 1)) x exp(program effects); library size
   decoupled from usage (row-normalized lambda). Replaces exp(Gaussian) lognormal
   counts. Usage noise 0.2 (A2) retained. Donor structure, truth, rare program, null
   columns, nested dropout, seeds unchanged.
2. **Gate operationalization (binding)**: because the pilot showed the absolute bars
   sit BELOW the linear ceiling for ANY method (SVD included), the registered gates
   are evaluated RELATIVE to the per-rep SVD ceiling, computed on the same normalized
   train matrix and stored per run:
   - A1 (operational): matched |cor| >= 0.9 x SVD-ceiling per planted program.
   - C2 (operational): rare-program |cor| >= 0.7 x SVD-ceiling in >= 8/10 seeds.
   The ABSOLUTE bars (0.9 / 0.7) are also reported for transparency, as is the
   method-vs-SVD efficiency gap (cebmf_cor - svd_ceiling). A method that matches the
   linear ceiling is doing everything a linear factorization can; a gate no method
   can meet measures nothing. The gates' INTENT (recover planted programs well) is
   unchanged; the operationalization makes it measurable on compositional data.
3. All other registered content (arms, metrics, protocol, seeds, gates C1/C3/C4,
   approach-B criteria) unchanged.
