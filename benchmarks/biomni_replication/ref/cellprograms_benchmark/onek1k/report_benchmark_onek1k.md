# OneK1K benchmark: cellprograms vs sample-level baselines

**Dataset**: OneK1K (healthy PBMC, 300-donor cap) → donor × cell-type pseudobulk, 9 cell types (43–300 donors per CT), top-3000 HVGs. Biology covariates: sex (categorical), Age (continuous). Technical: pool_number (73 pools / 300 donors). No cell-embedding comparators (GloScope/PILOT) for this cohort — embedding extraction exceeded memory and was not rerun.
**Method**: reduced ladder under the session time cap — L0 + L3 block on the L0 basis. Config: max_factors=25, seed 42.

## Headline results

1. **cp_L0 wins sex decisively: 0.978** — the strongest single Panel A result in the benchmark (grouped_pca 0.942, global_pca 0.849, composition 0.468, floor 0.518). Fourth cohort, fourth cp win or tie on sex; at 300 donors the kNN metric has the power to show the full separation.
2. **Age is compositional here**: composition 0.678 / clr 0.645 win; grouped_pca 0.556; cp_L0 0.402 (best per-program rep, above global_pca 0.336 and perct_pca 0.264). In a healthy cohort, aging's dominant blood signal is cell-composition shift (naive→memory), not within-cell-type programs — consistent with the compositional-phenotype pattern from COMBAT (Death28) and HLCA (smoking).
3. **Panel B: pools are well-mixed.** pool_number kNN near floor for all reps (cp_block 0.070, cp_L0 0.171, floor 0.006). The per-program screen flags 159/199 programs at η²>0.3, but with 73 pools over 300 donors the null expectation of η² is ≈(73−1)/(300−1) ≈ 0.24 — the screen is uninformative at this cardinality; trust the rep-level metric here.

## Panel A/B

![fig1](fig1_panelAB_heatmap.png)

## Panel C: representation quality and the L3 view

![fig2](fig2_r2_ladder.png)

L0 R² is 0.16–0.36 — markedly lower than the disease cohorts (COMBAT 0.35–0.85, HLCA 0.69–0.85, Stephenson 0.36–0.72). Two drivers: 300 donors × capped 25 factors (more independent biological variance to explain per factor), and a healthy cohort lacking the strong coordinated disease programs (e.g. interferon) that boost reconstruction R². Programs remain biologically real — witness the 0.978 sex separation.

![fig3](fig3_block_mass.png)

Block decomposition (L0 basis): 11 shared factors — **9 global, 2 partial** (NK-only; CD4_CM+CD4_Naive+CD8_EM+NK, a lineage-coherent lymphoid coalition).

## Panel E: practicality

| stage | minutes |
|---|---|
| fit_none (L0, 9 CTs, ≤300 donors) | 16.3 |
| block_none | 29.7 |
| **TOTAL** | **46.4** |

Fastest ladder of the four cohorts despite 2.5× the donors — donor count is not the runtime driver; L1 prior iterations are (fit_global was the dominant cost in COMBAT/HLCA and did not finish for Stephenson).

## Limitations

- Reduced ladder (L0 + block only): no L1/L2 comparison for this cohort.
- No GloScope/PILOT comparators (embedding not extracted).
- Low R² cohort — interpret programs via Panel A associations and cross-cohort matching rather than reconstruction alone.
- pool_number technical screen uninformative (high cardinality inflates η²).

## Artifacts

- `eval_all.csv`, `program_summary.csv`, `program_annotations.csv`, `block_structure.csv`, `runtime.csv` (shared `onek1k/fits/`)
- Figures: `fig1_panelAB_heatmap.png`, `fig2_r2_ladder.png`, `fig3_block_mass.png`
- Fits: `stage1_fits.rds` (L0), `integrations.rds` (block_none)
