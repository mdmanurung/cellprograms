# Cross-cohort immune-aging analysis with cellprograms

**Application report — Terekhova et al., Immunity 2026 (DOI 10.1016/j.immuni.2026.09.012)**
Run date: 2026-10-16. Scripts and committed tables: `analyses/immune_aging/`.
Case-study article: `vignettes/articles/immune_aging_cross_cohort.Rmd`.

## 1. Scope and design

The pipeline was applied to the 8-cohort PBMC immune-aging resource (2,609
healthy donors; 59 harmonized cell types): ABF300 (166 donors), OneK1K (981),
AIDA (625), CIMA (421), KCL (188), UCSF (99), SoundLife (96), SPAC (33).
Design (preregistered operating point, with one pilot-falsified amendment):

- **Pseudobulk**: sum of `.raw/X` integer counts over all cells per donor ×
  cell type (min 10 cells; "Bad quality" group dropped); log1p-CP10K.
- **Fits (MR-0)**: EBMF per cell type, point-Laplace prior, max_factors=10,
  top-2000 variable genes, var_type=1, seed=20261016, canonicalized.
  **No covariate residualization** — the original plan residualized
  age+sex+CMV inside the fit; a SPAC pilot showed residualized scores are
  orthogonal to age by construction (betas ~1e-5 SD/10y), making the
  association test vacuous. Amendment: fit unresidualized, test associations
  post hoc. Cohort-as-batch structure is handled by MR-0/MR-1 projection,
  not regression.
- **Associations**: `lm(score ~ age + sex + CMV)` per program (CMV dropped
  where constant), age beta per 10y in SD units, BH-FDR within cell type.
- **Projection (MR-1)**: discovery fit (ABF300; KCL for the CMV pair) frozen
  and projected onto each query cohort's pseudobulks with preprocessing
  replay; per-cell-type tryCatch isolation; min_coverage=0.9.
- **Conservation panel (3 tiers)**: (1) matched score concordance +
  association replication; (2) sign-flip/distribution stats; (3) top-20 gene
  Jaccard.

## 2. Program discovery and age prevalence

430–573 programs per cohort. Fraction of programs significantly
age-associated (FDR<0.05) scales with donor count:

| Cohort | Donors | Programs | Age FDR<0.05 | Fraction |
|---|---|---|---|---|
| CIMA | 421 | 573 | 320 | 0.558 |
| OneK1K | 981 | 571 | 249 | 0.436 |
| AIDA | 625 | 566 | 246 | 0.435 |
| UCSF | 99 | 430 | 83 | 0.193 |
| KCL | 188 | 564 | 106 | 0.188 |
| ABF300 | 166 | 561 | 48 | 0.086 |
| SoundLife | 96 | 552 | 52 | 0.094 |
| SPAC | 33 | 538 | 16 | 0.030 |

![Age prevalence](../vignettes/articles/figures/fig_aging_age_prevalence.png)

## 3. Thymic involution replicates across cohorts

Strongest-age program per cell type, cross-cohort sign consistency
(significant cohorts only):

| Cell type | Cohorts sig. | neg/pos | Median beta (SD/10y) |
|---|---|---|---|
| CD8+ RTE | 7 | 7/0 | −0.29 |
| CD4+ RTE | 7 | 6/1 | −0.32 |
| gd T naive | 7 | 6/1 | −0.32 |
| gd T Vd2 GZMK+ | 7 | 7/0 | −0.28 |
| CD8+ Temra | 7 | 5/2 | −0.22 |
| CD8+ Tem GZMK+ | 7 | 3/4 | +0.19 |
| CD8+ T Naive | 8 | 2/6 | +0.33 |
| CD4+ T Naive | 8 | 2/6 | +0.32 |

![RTE consistency](../vignettes/articles/figures/fig_aging_rte_consistency.png)

CD8+/CD4+ RTE program loadings are canonical recent-emigrant biology (SELL,
CD8B, TRBC2, CD8A, CD7; SOX4, LEF1). The naive-compartment *increase* is a
within-compartment shift (naive-vs-cytotoxic loadings), complementary to the
paper's proportion-level naive decline. The GZMK/GZMB health axis is a
proportion; at program level the Tem GZMK+ signal marks its
less-differentiated CD27+/LTB+ pole.

## 4. Cross-cohort transfer (tier 1: score concordance)

Matched score concordance (optimal assignment on |cor| between projected
discovery programs and query refit programs):

| Pair (query vs discovery) | Matched pairs | Median \|cor\| | frac ≥ 0.5 |
|---|---|---|---|
| SPAC (33) vs ABF300 | 186 | 0.589 | 0.602 |
| UCSF (99) vs ABF300 | 161 | 0.486 | 0.491 |
| AIDA (625) vs ABF300 | 184 | 0.591 | 0.576 |
| CIMA (421) vs ABF300 | 194 | 0.552 | 0.552 |
| OneK1K (981) vs ABF300 | 194 | 0.613 | 0.598 |
| SoundLife (96) vs KCL | 189 | 0.634 | 0.635 |

![Concordance](../vignettes/articles/figures/fig_aging_concordance.png)

Four rare KCL cell types fail projection (gd T Vd1 GZMB+/GZMK+, CD4+ Treg
cytotoxic, CD4+ Proliferative: rank-deficient shared-gene support) and are
logged, not fatal. OneK1K's "ILCs" type is absent from the ABF300 fit
(NA activities, warned).

## 5. Association replication (tier 1b)

**Matched-refit profiles are null-to-weak everywhere** (Spearman of age
betas across matched pairs: SPAC −0.07, UCSF 0.17, AIDA 0.02, CIMA 0.12,
OneK1K −0.02, SoundLife-vs-KCL age 0.13 / sex 0.23 / CMV −0.15). At 625–981
donors beta noise is negligible, so the failure is structural: unsupervised
EBMF decompositions distribute age variance across programs differently per
cohort (rotational freedom, adaptive shrinkage, greedy addition). Program
identity is cohort-specific even when content is conserved.

**Projected-score replication sidesteps the mismatch** — discovery programs
tested on their projected scores in each query:

| Query | Programs tested | Spearman (all) | Discovery hits replicated (p<0.05, same sign) |
|---|---|---|---|
| CIMA (421) | 558 | 0.25 | **68%** |
| OneK1K (981) | 551 | 0.14 | 52% |
| UCSF (99) | 528 | 0.12 | 41% |
| AIDA (625) | 548 | −0.10 | 40% |
| SPAC (33) | 534 | −0.10 | 36% |

Null ≈ 4%. AIDA is an unexplained outlier (projection quality, age range,
and ancestry comparable to CIMA; CI ±14% on 47 hits).

![Replication](../vignettes/articles/figures/fig_aging_replication.png)

## 6. Gene-set overlap (tier 3)

Top-20 gene Jaccard for matched program pairs: median 0.03–0.13, ≥0.2 in
8–36% of pairs. EBMF's dense, shrunk loadings make top-N gene identity
unstable even when the score subspace is conserved. Conservation claims
should be made at the score/subspace level.

![Jaccard](../vignettes/articles/figures/fig_aging_jaccard.png)

## 7. CMV axis (KCL → SoundLife)

The only CMV-informative pair. Matched score concordance 0.634; matched
CMV-profile replication null (Spearman −0.15, sign agreement 0.46). The
paper's CMV effect operates at the composition level (GZMK+/GZMB+ Tem
ratio); program-level CMV associations do not replicate under projection
matching. Composition-level analysis is future work.

## 8. Related registered benchmarks (same session)

- **Batch adjustment v1** (`benchmarks/batch_adjustment/RESULTS.md`): OLS
  residualization vs ComBat vs ComBat-seq vs fixed-covariate on simulated
  cohort-as-batch pseudobulks. Rule did not fire; OLS default retained.
- **Small-batch v1a**: 8 and 5 donors/batch. Rule did not fire in the gate
  scenario (ComBat leaves residual batchR2 0.028 vs OLS 0.000); S2 guidance
  documented (at 5 donors/batch OLS degrades: recovery 0.82, K inflated).
- **C2 matching check**: greedy vs optimal (Hungarian) program matching on
  the real K30 stability refits; pooled discordance 1.35% (Stephenson 1.19%,
  combat 1.58%) < 5% threshold → greedy matching retained, no R/ change.

## 9. Limitations

- Program sign orientation (canonicalization by top-|loading| gene) is
  mostly but not provably stable across fits; matched metrics use |cor|.
- CMV adjustment impossible outside KCL/SoundLife; ethnicity confounded
  with cohort (AIDA/CIMA/SPAC Asian; OneK1K/ABF300 European).
- SPAC (60+ only) has restricted age range; its age betas are noisy (n=33).
- Discovery limited to ABF300 (166 donors) — a larger discovery freeze
  (e.g. OneK1K) is a robustness check worth running.
- Tier-2 principal-angle subspace metrics and the stability-normalized
  conservation ratio were registered as extensions and are not yet computed.

## 10. Reproducibility

Environment: R 4.4.2, cellprograms v0.2.0+ (this repo), rhdf5, clue, sva;
Python 3.10+ (h5py, numpy, scipy). Data: open portal (access recipe in
`analyses/immune_aging/README.md`). Full command sequence, registered
operating point, and amendment rationale: `analyses/immune_aging/README.md`.
Large artifacts (fits, projection outputs) accompany the session results
archive (`aging_pbmc_repro/`).
