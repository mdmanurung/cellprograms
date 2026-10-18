# CD8T subset aging analysis: cEBMF programs vs age (program effect)

User directive: aging-associated gene programs WITHIN CD8 T cell subsets
(program effect, separated from composition). Registered pre-run (ledger
2026-10-17); recipe: author_celltype=='CD8 T cell', drop Bad quality, drop
zero-count; fine-label subsets; 99 donors, age 21-75 (median 33); CMV/BMI
constant ("No_info") -> model = score ~ standardized age + sex.

## Design

- 8 subsets kept (>=1,000 cells, >=20 donors): Naive (29,249 cells after QC),
  Tem GZMB+ (6,217), Tem GZMK+ (4,896), Tcm CCR4- , HLA-DR+ (98 donors),
  Temra (1,525), MAIT (1,312), RTE (94 donors).
- Two routes per subset: SEACells metacells (10k cap on Naive; 133/82/54/28/25/
  20/17/14 metacells) and donor pseudobulk; identical top-2000 HVG panel.
- cEBMF fits at the registered operating point (K=10, laplace/laplace,
  column-wise, backfitting, seed family 20261016).
- Age association: OLS beta (score ~ age_std + sex) + Spearman, BH-FDR over
  programs x subsets x routes (160 tests).

## Numerical incident resolved en route (documented deviation 7)

GZMKpos diverged on ALL 12 family seeds at every z-clip tier (10/5/3/2).
Root cause diagnosed (diag_gzmk.py / diag_gzmk2.py / diag_gzmk3.py): tau clamps
at the backend cap 1e12 (= 1/NUMERICAL_EPS) from iteration 0; with the
pipeline's 1e8 clamp the ELBO still improves until iteration ~21, then
pi0_F goes non-finite for ONE factor -> everything NaN. Mechanism:
single-donor-spike genes are exactly fittable at K=10 on 99 donors ->
residual variance -> 0 -> sebetahat -> 0 -> NaN inside the laplace EBNM
optimizer. Fix: sebetahat floored at 1e-4 before every laplace prior fit
(class-level guard, same layer as the registered tau cap). With the floor the
fit converges (ELBO plateau, k_eff=10); well-conditioned coordinates
(se >> 1e-4) are unaffected. The z-clip escalation ladder (deviation 6) is
retained but was NOT the operative remedy for this subset.

## Results (all 8 subsets x 2 routes converged; fits.npz + CSVs)

19 FDR<0.05 associations of 160 tests (`age_associations.csv`; heatmap
`age_beta_heatmap.png`, grey = pruned programs). Highlights with top-8 loadings:

- **Naive / pseudobulk P1, beta=-6.88, FDR=6.9e-5 (strongest)**: CST7, MYO1F,
  KLRB1, PRF1, ARID5B, PLEK, ERN1, JUNB - an effector/cytotoxic-like program
  (granzyme K/perforin axis) DECLINING with age within naive CD8T.
- **Tem GZMB+ / pseudobulk P3, beta=+2.14, FDR=0.002**: GZMA, GZMH, ANXA1,
  SH3BGRL3 - cytotoxic granzyme-A/H program INCREASING with age in GZMB+ Tem.
- **Tem GZMK+ / pseudobulk P3, beta=-2.09, FDR=0.001**: ribosomal block
  (RPL5/RPS8/MT-ND3...) - see caveat below.
- **RTE / metacell P10, beta=+0.92, FDR=0.018**: CFH, CTSV, ROR1, DISP3 -
  ROR1/CTSV ontogeny-associated program increasing with age (recent thymic
  emigrants in older donors).
- Tcm CCR4-: P1 down / P2 up in BOTH routes (sign-consistent); P3
  (NIBAN1, NR4A2 - activation/tolerance) down.
- HLA-DR+ / metacell P3/P4/P10 all declining with age.

**Interpretation caveat (honest)**: several FDR-significant programs are
dominated by ribosomal/mitochondrial/housekeeping blocks (RPL/RPS/MT-, B2M,
TMSB4X) - likely content/composition axes rather than immune biology. The
biologically interpretable associations are the Naive effector decline, the
GZMB+ Tem cytotoxic increase, and the RTE ROR1/CTSV program. Program-level
age effects coexist with composition effects (`composition.csv`).

## Route concordance and sharing

- Route concordance (`route_concordance.csv`): mean matched |cor| 0.10-0.17;
  the FIRST factor matches strongly (0.60-0.77 in 5/8 subsets), remaining
  factors are route-specific - expected at K=10 over-parameterization.
- Sharing sensitivity (`sharing_sensitivity.csv`, real_data_case_study
  pattern): all 28 subset pairs pass the permutation calibration at both
  n_perm=199 and 1999 (28/28, ARI=1.0) - the null is tight (same
  loose-calibration behavior documented in workstream E). But at the
  |cor|>=0.8 magnitude bar there are ZERO shared edges on the 884-gene
  universe; only 2 (metacell) / 5 (pseudobulk) edges at >=0.7. Conclusion:
  cross-subset program sharing is statistically detectable but weak in
  magnitude - CD8T subset programs are largely distinct; the 199-vs-1999
  calibration is stable but uninformative at this effect size.

## Registered deviations (also in _run_stats.json)

1. 10k-cell seeded subsample for Naive metacell build; 2. zero-variance gene
drop (union across routes); 3. tau precision capped at 1e8; 4. seed-family
retries (first non-diverging seed wins); 5. z-clip [-10,10]; 6. z-clip
escalation ladder [-5,5]->[-3,3]->[-2,2] (both routes refit together);
7. sebetahat floored at 1e-4 before laplace EBNM fits (GZMKpos mechanism,
above). Age+sex-only model (CMV/BMI constant).

## Files

age_associations.csv, route_concordance.csv, sharing_sensitivity.csv,
composition.csv, fits.npz, age_beta_heatmap.png, _run_stats.json.
