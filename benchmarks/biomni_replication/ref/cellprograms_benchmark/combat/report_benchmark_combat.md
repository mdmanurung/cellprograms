# COMBAT benchmark: cellprograms vs sample-level baselines and comparators

**Dataset**: COMBAT COVID-19 blood atlas (CELLxGENE Census 2025-11-08, dataset ebc2e1ff) — 124 donors, 10 cell types, donor×cell-type pseudobulk (≥50 cells/sample), top 3000 HVGs per cell type, logCPM.
**Cohort composition**: COVID-19 severe (30), sepsis (23), COVID critical (17), mild (17), HCW-mild (13), flu (12), healthy (10), COVID-LDN (2). Technical covariates: Institute (Oxford 110 / St George's 14), GEX region.
**Pipeline under test**: per-cell-type empirical-Bayes matrix factorization (flashier) with the local-first borrowing ladder — L0 (no borrowing), L1-global (shared prior), L1-groups (lineage-pooled priors), L2 (initialization borrowing), L3 block joint decomposition, plus EV-BIDIFAC-S and PVA as integration alternatives.
**Comparators**: random embedding, cell-type composition, CLR composition, global pseudobulk PCA, per-CT PCA, grouped (lineage) pseudobulk PCA, GloScope (GMM density, KL), PILOT (Wasserstein).

## Headline results

1. **cellprograms representations are the cleanest biology carriers on unambiguous expression-driven traits** — sex (cp_L2 0.86 vs global_pca 0.82, composition 0.52) and Age (cp_L1groups 0.46 vs per-CT PCA 0.42, global PCA 0.30) — while sitting **at the random floor on the technical covariate** (Institute 0.47 = floor for all cp variants; per-CT PCA 0.67, grouped PCA 0.60).
2. **Cheap PCA baselines win the coarse disease axis** (Source: grouped 0.46, per-CT 0.45 vs cp ~0.28). COMBAT's disease signal is strongly compositional and global-interferon-driven; kNN in a 266-dimensional sparse program space dilutes it relative to a 10–30-dimensional dense PCA space. This is a metric–dimensionality interaction as much as a biology difference (see Limitations).
3. **L1-global ≈ L1-groups on Panel A** (max Δ = 0.039 on sex; identical on Death28/Hospitalstay/TimeSinceOnset). Program-level: 251/281 programs match across the two fits (median |cosine| 0.967); divergence concentrates in low-coverage CTs (GDT, PB, DP) on continuous clinical covariates, bidirectionally.
4. **L2 initialization borrowing damages the rare cell types it targets**: PB 25→5 programs (R² 0.57→0.13), GDT 24→8 (R² 0.56→0.20), DP 18→12 (R² 0.35→0.25). A genuine negative result from the real-data test.
5. **L3 block decomposition recovers lineage-coherent sharing**: 9 global + 11 partial factors; partials respect lineage boundaries (CD8+DP, cMono+ncMono+DC, CD4+NK, CD4+CD8+NK+cMono ×3). Identical classification on L0 and L1-global bases.
6. **EV-BIDIFAC-S is technically contaminated** (Institute 0.87 — the worst of any representation) and its apparent Panel A wins (Source 0.50) are confounded by batch.

## Panel A/B: covariate preservation

![Panel A/B heatmap](fig1_panelAB_heatmap.png)

Values: leave-one-out 5-NN macro-F1 (categorical) / max |Spearman| across dims (continuous). Random floor included as first row.

Reading:
- **sex** is the cleanest win for the program representations (0.77–0.86; all four ladder variants ≥ global PCA). Sex programs are recovered as discrete sparse programs (Y-chromosome + XIST; see Panel D) in every cell type.
- **Age** favors cp variants (0.43–0.46) over all baselines (≤0.42) and both comparators (≤0.16).
- **Source / Outcome / Death28** favor composition-aware baselines. Death28 is near-purely compositional (composition 0.57; all expression-based reps 0.46–0.52).
- **PILOT underperforms throughout** (≤0.23 on continuous covariates); **GloScope-GMM** is mid-pack with a compositional tilt (Death28 0.57).
- **Panel B**: all cp ladder variants and cp_block sit exactly at the Institute floor (0.47); cp_evbS is an outlier at 0.87 — its joint states absorb the site batch, which also inflates its Source/Outcome scores (Institute and Source are confounded: St George's processed specific clinical groups). GEX_region is near-floor for everything.

## L1-global vs L1-groups: where they diverge

![L1 divergence](fig4_l1_divergence.png)

- Representation level: max |Δ| = 0.039 (sex). The two donor-level representations are effectively equivalent on Panel A.
- Program level: 251 matched pairs, median |cosine| 0.967; 15 programs unique to each side (borderline factors swapped at the 30-factor cap in major CTs).
- Divergence concentrates in **low-coverage CTs** (GDT, PB, DP, B) and on **continuous clinical covariates** (Age, Hospitalstay, TimeSinceOnset; max |Δ| ≈ 0.27). Direction is bidirectional — lineage pooling recovers signal in some programs (GDT_F009 Age: 0.06→0.32) and attenuates it in others (B_F024 TimeSinceOnset: 0.25→0.01).
- **Interpretation**: on a cohort of this size the pooling-scope choice is immaterial at the representation level and acts only at the sparse-data margins. L1-groups remains the recommended default on principle (it cannot homogenize across lineages — the failure mode L1-global showed in simulation), not on a measured Panel A advantage here.

## Panel C: representation quality and the L3 view

![R2 ladder](fig2_r2_ladder.png)

- Per-CT R² is stable across L0/L1-global/L1-groups (±0.01); L1 variants slightly *reduce* program counts in rare CTs (PB 25→23, GDT 24→22, DP 18→17) at equal R² — the intended regularization.
- **L2 collapses on its target CTs** (PB, GDT, DP): initialization from the donor CT's programs traps the fit in a poor local optimum. L2 as currently implemented should not be used for rare-CT rescue; the ladder's L1 rungs are the safe borrowing mechanism. (Its donor-level Panel A scores remain competitive because the shared major-CT programs dominate the wide representation.)

![Block mass](fig3_block_mass.png)

- Block decomposition (expression level, 20 shared factors): 9 global, 11 partial. Partials are lineage-coherent: CD8+DP (SF06/07; DP are double-positive T cells), myeloid axis (SF08 cMono+ncMono+DC, SF10 ncMono+DC, SF20 ncMono), and a repeated CD4+CD8+NK+cMono activation coalition (SF11/16/19 — the interferon/inflammation space; see Panel D).
- **SF18 (DC+PB+GDT+DP)** spans the four lowest-coverage CTs and carries the largest single-factor PVE (0.167). Flagged for caution: plausibly plasmablast proliferation amplitude, possibly low-coverage noise alignment.
- The classification is **identical on L0 and L1-global bases** — prior pooling did not shift any factor's sharing class.
- Private residual pass saturated at 10 factors/CT (default cap): substantial CT-private structure remains beyond the 20 shared factors.

## Panel D: interpretability

Program annotations (L1-groups fit, 266 programs; top-10 genes by |loading|, Ensembl→symbol via dataset feature metadata):

- **Sex programs** recur in every CT with the expected gene set (XIST, RPS4Y1, DDX3Y, UTY, KDM5D, EIF1AY): e.g. CD4__F008, CD8__F009, B__F006.
- **Type-I interferon programs** recur across CTs (IFIT3, IFI44L, RSAD2, IFIT1, USP18, OAS1): e.g. CD4__F006, CD8__F010, B__F005, PB__F003; the monocyte instance (cMono__F002) is led by SIGLEC1, a known IFN-induced monocyte marker in COVID-19. These cross-CT recurrences are the program-level counterpart of the block decomposition's global factors.
- **Plasmablast dominant program** (PB__F001, PVE 0.11) is a proliferation signature (ESPL1, CDC25A, CEP15, TRAIP) — consistent with plasmablast blast biology.
- **Loading concentration**: median top-50 |loading| mass = 0.10 — programs are distributed pathway-scale objects, not single-gene modules.
- **Technical contamination screen**: 2/266 programs with Institute η² > 0.3 (CD4__F007 0.31, cMono__F004 0.37); 0/266 for GEX_region. The L1-groups fit is essentially clean at program level, consistent with the floor-level Institute kNN score.
- **Sharing spectra (L4)**: loadings space (gene usage, residual-SD whitened) shows the expected myeloid coherence (cMono↔ncMono: min angle 10°, 14 shared directions) and broad low-level sharing across pairs (45/45 pairs ≥1 direction at ~1000-dim effective gene space). **Scores-space spectra are uninformative on this dataset**: with ~30 programs per CT and only ~50 shared donors per pair, subspace dimensions exceed the donor overlap (k_a + k_b > n), so principal angles degenerate to 0 by construction. Use loadings space at this program density; scores-space sharing requires either fewer programs per CT or more shared donors.

## Panel E: practicality

| Stage | Minutes |
|---|---|
| fit_none (L0) | 20.1 |
| fit_global (L1) | 364.4 |
| fit_groups (L1) | 98.0 |
| fit_l2 (L2) | 19.7 |
| integrations (block ×2, PVA, EV-BIDIFAC) | 75.7 |
| **Total ladder** | **578.9** |

- L1 prior pooling is the dominant cost (L1-global ≈ 18× L0 at prior_max_iter=5; two CTs hit flashier's iteration ceiling). Rollout datasets use the trimmed config (max_factors=25, prior_max_iter=3).
- Comparator practicality notes: GloScope's KNN density estimator produced pervasive non-finite distances on this data (117/124 donors affected, KL and JS alike — zero kNN density → log 0); the GMM variant is stable (finite, ~2 min). PILOT is fast (~10 s) but weak on every covariate here.

## Stability

Subsample stability (10 refits at 80% donor resampling, L0 fit, 6 major CTs, 180 programs, correlation matching):

- **132/180 programs (73%) form a stable core** — recovered in ≥8/10 subsamples — with median loading cosine 0.854 against the full-data program.
- Among recovered programs: median score correlation 0.855, median loading cosine 0.818. Top-gene Jaccard is modest (median 0.212) — expected for distributed programs (median top-50 mass 0.10): the subspace is stable even when the exact top-gene ranking varies.
- Stability tracks program importance: PVE vs recovery frequency Spearman ρ = 0.40. The 28 never-recovered programs are the low-PVE tail — minor factors that come and go, not the programs one would interpret.
- Per cell type: CD4/CD8/cMono most stable (83–87% stable core), B and ncMono least (53–60%).

Interpretation: the interpretable program set is robust to donor resampling; instability is confined to low-variance minor factors.

## Limitations

1. **Metric–dimensionality interaction**: kNN macro-F1 in 231–273-dimensional program space penalizes cp representations relative to 10–30-dimensional PCA embeddings on coarse axes. Per-program association tests (Panel D, divergence table) are the dimensionality-fair complement and show cp programs do carry disease-associated signal.
2. **Composition–expression confounding**: COMBAT's severity/mortality axes are strongly compositional; no expression representation was designed to win those, and composition should be reported alongside programs in practice.
3. **Missingness**: only 9/124 donors have all 10 CTs at ≥50 cells; donor-aligned representations use 0-imputation after z-scoring (documented in methods).
4. **Scores-space L4 degeneracy** at high program counts (above).
5. Single-cohort technical structure (2 sites, imbalanced 110/14) limits Panel B resolution.

## Artifacts

- `fits/`: stage1_fits.rds, integrations.rds, block_none_preview.rds (L0 + block cross-check)
- Tables: eval_all.csv, l1_divergence_programs.csv, l1_divergence_by_ct.csv, program_summary.csv, program_annotations.csv, block_structure.csv, sharing_spectrum_scores.csv, sharing_spectrum_loadings.csv, runtime.csv, pilot_dist.csv, gloscope_gmm_dist.rds
- Figures: fig1_panelAB_heatmap.png, fig2_r2_ladder.png, fig3_block_mass.png, fig4_l1_divergence.png
