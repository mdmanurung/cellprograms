# PREREG batch-adj-v1: batch-effect adjustment for pseudobulk EBMF programs

Status: REGISTERED (draft committed before any benchmark run; results land in
RESULTS.md exactly as specified here — results recorded either way, adoption
only if the rule fires).

## Question

For donor-level pseudobulk EBMF (cellprograms), does empirical-Bayes batch
adjustment (sva::ComBat on log-pseudobulk; sva::ComBat-seq on count-pseudobulk)
or within-model batch factors (covariate_mode = "fixed") recover planted
programs better and remove more batch variance than the current default
(OLS residualization, covariate_mode = "residualize")?

## Simulator

Base: `simulate_pseudobulk(scenario = "mixed", n_samples = 150, n_genes = 1500,
program_size = 50, signal_sd = 1.5, noise_sd = 1.0)` — 5 cell types, planted
programs with known donor activities z (truth).

Batch planting (extension, benchmarks-only code): assign each donor to one of
B = 4 batches (balanced). Per batch, shift a random 10% gene subset by
delta_g ~ N(0, sigma_b) on the log scale, sigma_b in {0.375, 0.75} (0.25x and
0.5x signal_sd), plus a per-batch scale multiplier on all genes
(1 +/- 0.05). Batch labels are independent of program activities z (batch is
not confounded with truth in the main arm; a confounded arm — 30% |cor(batch,
z)| — is included as a secondary setting because real cohorts have it).

Count arm: exponentiate the adjusted log-profiles with donor library sizes
(log-normal, sigma 0.4), round to integer counts, run ComBat-seq
(batch = batch), re-log-normalize (library-size + log1p) before fitting.

## Arms (matched design; one factor changes per comparison)

- A0 (default): OLS residualization of batch dummies —
  `fit_celltype_programs(covariates = ~batch)` (residualize).
- A1: `sva::ComBat` on log-pseudobulk (batch = batch, no preserved covariates),
  then fit with no covariates.
- A2: `sva::ComBat-seq` on count-pseudobulk, re-normalized, fit with no
  covariates (counts arm only).
- A3: within-EBMF fixed factors — `covariate_mode = "fixed"`,
  `covariates = ~batch`.

All fits: max_factors = 10, point_laplace, var_type = 1, features = "all"
(simulator plants programs across the full gene panel), seed = 1, canonicalize
after fit. 10 simulation seeds per sigma_b per arm; cell type "B" and "CD8"
analyzed (one private-program type, one shared-program type).

## Metrics (primary, per seed x sigma_b x cell type)

1. Program recovery: mean over planted programs of max |cor(canonical donor
   scores, truth z)| (greedy K may return extra factors; take the max over
   factors, then average over programs).
2. Batch variance removed: 1 - R^2 of each recovered program's scores on batch
   dummies, averaged (higher = less batch in scores).
3. Secondary: stability (existing refit correlation, 5 refits, one setting);
   CMP membership F1 at the CMP-detection operating point (one setting).

## Decision rule (registered)

Adopt an EB mode (A1/A2) or the fixed mode (A3) into R/ ONLY IF, pooled over
seeds and sigma_b at the primary setting: recovery improves by >= 0.02 mean
correlation vs A0 (paired per seed) AND batch variance in scores is >= 10%
lower than A0, with no stability metric degradation > 0.02. Otherwise the
default stays and the result is recorded. Confounded-arm results are reported
but do not gate adoption (they inform the real-cohort interpretation).

## Real-data sanity (non-gating)

One cohort with CMV metadata (KCL or Sound Life, per the Phase 1 audit): fit
all feasible arms; compare program-score associations with age and CMV status
(direction and magnitude). Reported in RESULTS.md as context.

## Power/precision note

10 seeds x 2 sigma_b x 2 cell types x 4 arms = 320 fits at ~seconds each
(150 donors x 1500 genes): minutes of compute. Paired per-seed comparison
detects a 0.02 correlation shift at ~0.8 power given the per-seed SD observed
in prior simulator benchmarks (~0.03); if the observed SD exceeds 0.06, extend
to 20 seeds before reading the rule (pre-registered extension).
