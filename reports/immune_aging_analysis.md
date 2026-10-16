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

### 4a. Discovery-cohort sensitivity (OneK1K refreeze)

The whole panel was recomputed with **OneK1K (981 donors) as the discovery
freeze** (MR-0 refit at the registered operating point; MR-1 projection of
all 7 remaining cohorts). Per-cohort comparison with the registered
ABF300-discovery panel (`tables/sensitivity/discovery_sensitivity.csv`):

| Query | Median \|cor\|, ABF300 disc. | Median \|cor\|, OneK1K disc. | Sign flips (ABF3 / OneK1K) |
|---|---|---|---|
| ABF300 (166) | — | 0.545 | — / 0.330 |
| SPAC (33) | 0.589 | 0.528 | 0.392 / 0.363 |
| UCSF (99) | 0.486 | 0.579 | 0.304 / 0.343 |
| AIDA (625) | 0.591 | 0.538 | 0.321 / 0.315 |
| CIMA (421) | 0.552 | 0.634 | 0.361 / 0.329 |
| KCL (188) | — | 0.508 | — / 0.260 |
| SoundLife (96) | 0.634 (vs KCL) | 0.710 | 0.365 / 0.291 |

**Verdict: the conservation panel is discovery-robust.** Switching the
discovery freeze moves matched score concordance by at most ±0.09 per
cohort (OneK1K-disc range 0.508–0.710 vs ABF300-disc 0.486–0.634), slightly
widening the panel's upper end (SoundLife 0.710) without collapsing any
pair. Matched age-profile Spearman stays null-to-weak under both
discoveries (OneK1K-disc: −0.079..0.116; ABF300-disc: −0.07..0.17) — the
decomposition-mismatch conclusion is a property of the method, not of the
discovery cohort.

Projected-score age replication under the OneK1K discovery
(`tables/sensitivity/projected_age_replication_onek1k_disc.csv`): the
981-donor freeze yields 244–249 significant age hits (FDR<0.05; ABF300's
166-donor freeze yielded 45–48), and 16–37% of them replicate in each query
cohort (p<0.05, same sign; ABF300 0.198, SPAC 0.262, UCSF 0.219, AIDA 0.269,
CIMA 0.365, KCL 0.325, SoundLife 0.162; null ~0.04). Matched age-profile
Spearman spans −0.03..0.46 (highest KCL 0.462, SPAC 0.314). The larger
discovery's hit set includes smaller effects that a <200-donor query is
underpowered to confirm, so the lower raw replication fraction vs the
ABF300-discovery numbers (36–68% of 45–48 hits) reflects hit-set depth, not
weaker transfer — the qualitative conclusion (partial score-level transfer,
well above null) is unchanged.

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

## 7. Composition-level aging axis (proportions and the Tem ratio)

Program scores capture within-compartment transcriptional shifts; cell
**proportions** capture a second, independent aging axis. Per-donor
per-cell-type counts were extracted from each cohort h5ad and modeled
directly (`prop ~ age + sex [+ CMV]` per key type; scripts
`composition_counts.py` + `composition_aging.R`, tables under
`tables/composition/`).

**Proportion-level thymic involution is universal** — more consistent
across cohorts than the program-level declines (Section 3):

| Key type | Cohorts significant (p<0.05) | Relative change per decade |
|---|---|---|
| CD8+ RTE | 8/8 declining | −22% to −70% |
| CD4+ RTE | 8/8 declining | −20% to −63% |
| CD8+ T Naive | 8/8 declining | −19% to −42% |
| gd T Vd2 GZMK+ | 7/8 declining | −8% to −29% |
| CD8+ Tem GZMB+ | 6/8 increasing | +9% to +18% |
| CD8+ Tem GZMK+ | 6/8 increasing | +4% to +23% |
| CD8+ Temra | 3/8 increasing | +9% to +15% |
| CD4+ T Naive | 3/8 (mixed) | ±3% |

![Composition heatmap](../vignettes/articles/figures/fig_aging_composition.png)

**The GZMK+/GZMB+ Tem ratio is not a monotonic age effect.** Its age
trajectory is cohort-inconsistent — log-ratio slope per decade −0.101
(OneK1K, p = 7.6e-12), −0.078 (CIMA, p = 1.6e-3), −0.161 (SPAC, p = 1.2e-3),
but **+0.157 in AIDA (p = 5.4e-6, opposite sign)** and null in ABF300, KCL,
UCSF, SoundLife. This matches the paper's framing of the ratio as a
health-trajectory stratifier rather than an age effect.

| Cohort | n | Median ratio | Age slope (log/10y) | age p | CMV effect (log) | CMV p |
|---|---|---|---|---|---|---|
| ABF300 | 166 | 1.041 | +0.077 | 0.16 | — | — |
| OneK1K | 981 | 0.590 | −0.101 | 7.6e-12 | — | — |
| AIDA | 625 | 0.622 | +0.157 | 5.4e-6 | — | — |
| CIMA | 421 | 0.399 | −0.078 | 1.6e-3 | — | — |
| KCL | 188 | 0.761 | +0.034 | 0.29 | −0.31 | 0.61 |
| UCSF | 99 | 0.856 | −0.043 | 0.40 | — | — |
| SoundLife | 96 | 0.991 | −0.012 | 0.80 | **−1.07** | **2.5e-11** |
| SPAC | 33 | 0.599 | −0.161 | 1.2e-3 | — | — |

**The CMV effect on the ratio replicates cleanly** where serostatus is
recorded: in SoundLife, CMV+ donors shift the log ratio by −1.07
(p = 2.5e-11) — a ~66% lower GZMK+/GZMB+ ratio in CMV+ donors, exactly the
paper's claim. KCL's smaller, null estimate (−0.31, p = 0.61) is diluted by
its No_info stratum (kept as its own level).

![Tem ratio](../vignettes/articles/figures/fig_aging_tem_ratio.png)

## 8. CMV axis at program level (KCL → SoundLife)

The only CMV-informative pair. Matched score concordance 0.634; matched
CMV-profile replication null (Spearman −0.15, sign agreement 0.46). The
paper's CMV effect operates at the composition level (GZMK+/GZMB+ Tem
ratio, Section 7); program-level CMV associations do not replicate under
projection matching.

## 9. Related registered benchmarks (same session)

- **Batch adjustment v1** (`benchmarks/batch_adjustment/RESULTS.md`): OLS
  residualization vs ComBat vs ComBat-seq vs fixed-covariate on simulated
  cohort-as-batch pseudobulks. Rule did not fire; OLS default retained.
- **Small-batch v1a**: 8 and 5 donors/batch. Rule did not fire in the gate
  scenario (ComBat leaves residual batchR2 0.028 vs OLS 0.000); S2 guidance
  documented (at 5 donors/batch OLS degrades: recovery 0.82, K inflated).
- **C2 matching check**: greedy vs optimal (Hungarian) program matching on
  the real K30 stability refits; pooled discordance 1.35% (Stephenson 1.19%,
  combat 1.58%) < 5% threshold → greedy matching retained, no R/ change.

## 10. Limitations

- Program sign orientation (canonicalization by top-|loading| gene) is
  mostly but not provably stable across fits; matched metrics use |cor|.
- CMV adjustment impossible outside KCL/SoundLife; ethnicity confounded
  with cohort (AIDA/CIMA/SPAC Asian; OneK1K/ABF300 European).
- SPAC (60+ only) has restricted age range; its age betas are noisy (n=33).
- Discovery was registered at ABF300 (166 donors); the OneK1K (981-donor)
  refreeze (Section 4a) confirms the panel is discovery-robust, but only
  two of the eight cohorts have served as discovery so far.
- Tier-2 principal-angle subspace metrics and the stability-normalized
  conservation ratio were registered as extensions and are not yet computed.

## 11. Reproducibility

Environment: R 4.4.2, cellprograms v0.2.0+ (this repo), rhdf5, clue, sva;
Python 3.10+ (h5py, numpy, scipy). Data: open portal (access recipe in
`analyses/immune_aging/README.md`). Full command sequence, registered
operating point, and amendment rationale: `analyses/immune_aging/README.md`.
Large artifacts (fits, projection outputs) accompany the session results
archive (`aging_pbmc_repro/`).
