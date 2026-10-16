# RESULTS: batch-effect adjustment benchmark (batch-adj-v1)

Registered in `PREREG.md` (committed before any benchmark run). Simulator:
`simulate_pseudobulk` (150 donors, 1500 genes, 2 planted programs of size 50,
signal_sd 1.5, noise_sd 1.0) + planted batch: 4 balanced batches, 10% of genes
shifted by delta ~ N(0, sigma_b), per-batch global scale shift +/-0.05; a
confounded variant matches batch quantiles on the first planted program
(a = 0.3, non-gating). Arms: **A0** OLS residualization on batch dummies
(package default, `covariate_mode = "residualize"`); **A1** `sva::ComBat` on
log-pseudobulk, no covariates after; **A2** `sva::ComBat-seq` on count
pseudobulk, no covariates after; **A3** `covariate_mode = "fixed"`
(within-EBMF fixed batch factors). 10 seeds x 2 cell types (B, CD8) per
scenario. Metrics: planted-program recovery (mean max |cor(score, truth z)|),
residual batch R^2 of donor scores on batch dummies (lower is better), number
of programs K.

## Results

| Scenario | Arm | Recovery | Batch R^2 | K |
|---|---|---|---|---|
| sigma_b = 0.375 | A0 OLS | 0.985 | ~0 (2e-9) | 1.5 |
| | A1 ComBat | 0.992 | 4.6e-3 | 1.5 |
| | A2 ComBat-seq | 0.976 | 1.8e-2 | 2.5 |
| | A3 fixed | 0.985 | ~0 (5e-7) | 1.5 |
| sigma_b = 0.75 | A0 OLS | 0.985 | ~0 (2e-9) | 1.5 |
| | A1 ComBat | 0.989 | 2.0e-3 | 1.5 |
| | A2 ComBat-seq | 0.976 | 4.8e-1 | 5.5 |
| | A3 fixed | 0.985 | ~0 (5e-7) | 1.5 |
| confound a = 0.3 (non-gating) | A0 OLS | 0.980 | ~0 (3e-9) | 1.5 |
| | A1 ComBat | 0.990 | 6.3e-3 | 1.5 |
| | A2 ComBat-seq | 0.976 | 2.1e-2 | 2.5 |
| | A3 fixed | 0.980 | ~0 (1e-6) | 1.5 |

Paired per-seed deltas vs A0 (sigma = 0.375): A1 recovery +0.0068 (sd 0.0043),
batch R^2 +0.0046 (i.e. ComBat leaves residual batch that OLS removed); A2
recovery -0.0086 (sd 0.0061), batch R^2 +0.0179, K inflated; A3 recovery
+0.0000 (sd 0.0001), batch R^2 +0.0000 (bit-equivalent to A0 when batch is
orthogonal to truth). Same pattern at sigma = 0.75 (A2 batch R^2 0.48, K 5.5)
and under confounding (A1 recovery +0.0101, still leaves batch R^2 6.3e-3).

## Registered decision

Adopt an EB mode into `R/` only if recovery improves by >= 0.02 (paired mean)
AND batch R^2 is >= 10% lower than A0, without stability degradation.

**The rule does not fire; the package default (OLS residualization) is
retained.** Neither arm meets either primary condition in any scenario:

- A1 (ComBat): recovery gain +0.007 to +0.010 (needs +0.02); batch R^2 is
  *higher* than A0 (EB shrinkage under-corrects at these batch sizes).
- A2 (ComBat-seq): recovery is *lower* (-0.009), batch R^2 is higher, and K is
  inflated (batch variance re-emerges as extra programs).
- A3 (fixed): numerically identical to A0 when batch is orthogonal to the
  planted programs; no gain, as predicted analytically.

Per the registered extension rule (per-seed SD > 0.06 triggers 20 seeds): all
per-seed SDs were 0.003-0.011, so no extension was needed.

## Secondary metric: stability

The harness's stability extraction was buggy (it unlisted the fit object
instead of reading `$stability$per_program`), so stability values are NA in the
saved results. The harness code is fixed for future runs. The decision was
decided by the primary criteria, which fail decisively in both directions
(recovery gain too small AND residual batch higher), so the missing secondary
metric cannot change the outcome. No `R/` change is made; no NEWS entry
required per repo law.

## Interpretation

With batch orthogonal to biological truth and balanced designs, OLS
residualization removes batch exactly (it is the correct generative model),
while EB shrinkage deliberately shrinks correction parameters and therefore
leaves residual batch in the data. ComBat's advantage appears with small or
unbalanced batches; ComBat-seq's with count-level overdispersion. Neither
condition holds in the target use case (donor-level pseudobulks, balanced
cohort-as-batch designs handled by MR-0/MR-1 projection rather than
regression). The confounded scenario (a = 0.3) shows the same ordering and is
non-gating by registration.

Reproducibility: `run_benchmark.R --seeds 10 --sigma 0.375|0.75 --confound 0|0.3`;
raw results in `results.rds` per scenario directory.
