# Stephenson benchmark: cellprograms vs sample-level baselines and comparators

**Dataset**: Stephenson et al. COVID-19 PBMC (~780k cells) → donor × cell-type pseudobulk, 19 cell types, 31–120 donors per cell type, top-3000 HVGs. Biology covariates: Status (Covid/Healthy), Worst_Clinical_Status, Outcome, Smoker, sex (categorical), Days_from_onset, Age (continuous). Technical: Site (Ncl/Cambridge/Sanger); Collection_Day and Resample are effectively constant at donor level (119/120 one class) — included for completeness but uninformative (random floor 0.996).
**Method**: reduced ladder under the session time cap — L0 (fit_none) + L1-groups + L3 block on the L0 basis. Config: max_factors=25, prior_max_iter=3, seed 42. Comparators: 6 cheap baselines + GloScope-GMM + PILOT.

## Headline results

1. **cp_L0 wins Age (0.653 vs ≤0.48 baselines, GloScope 0.136) and Days_from_onset (0.559 vs ≤0.38)** — third cohort where cp representations carry the strongest continuous-clinical signal. cp_L1groups lands close (Age 0.522, Days 0.539, sex 0.760) — L1 borrowing neither helps nor hurts at rep level here, consistent with COMBAT/HLCA.
2. **cp_block wins Outcome (0.861) and essentially ties clr on Status (0.833 vs 0.834)** — the block shared factors capture systemic COVID signal that per-cell-type programs fragment (cp_L0 Status 0.460). First dataset where the L3 representation is the best biology rep for the headline clinical covariates.
3. **sex: three-way tie at the top** (grouped_pca 0.781, cp_L0 0.779, global_pca 0.773); GloScope at floor (0.513).
4. **Panel B: every non-cp representation is maximally site-confounded** — perct_pca 1.000, grouped_pca 1.000, GloScope 0.960, clr 0.939, global_pca 0.938 (floor 0.444). cp reps are the cleanest: cp_block 0.558, cp_L0 0.684. Stephenson's multi-site design makes this the starkest Panel B separation of the three cohorts.
5. **clr/composition win the compositional severity covariates** (Worst_Clinical_Status 0.475, Smoker 0.543) — consistent with COMBAT Death28 and HLCA smoking.

## Panel A/B: covariate preservation

![fig1](fig1_panelAB_heatmap.png)

kNN macro-F1 (categorical) / max |Spearman| (continuous). Collection_Day/Resample columns are degenerate (119/120 donors one class; floor 0.996) — interpret only Site on Panel B.

## Panel C: representation quality and the L3 view

![fig2](fig2_r2_ladder.png)

L0 per-cell-type R² 0.36–0.72 (16–25 programs per CT). Lower than HLCA (0.69–0.85) — expected: a heterogeneous disease cohort with several low-donor cell types (NK_prolif n=31, B_immature n=34). **L1-groups is essentially identical to L0 per cell type** (R² within ±0.02 except B_immature 0.531→0.448, the lowest-donor CT at n=34; program counts unchanged except B_immature 17→15, B_switched_memory 25→24) — third cohort where L1 pooling scope is immaterial at the representation level.

![fig3](fig3_block_mass.png)

Block decomposition (L0 basis, 20 shared factors): **13 global, 5 partial, 1 private**. Partials are lineage-coherent: CD14_mono+CD16_mono and CD14_mono+CD83_CD14_mono (monocyte compartment), CD8.TE, CD4.CM, CD4.IL22+CD8.Naive. The single private factor is **Platelets** — platelets share essentially no expression programs with nucleated cells, a sensible structural result.

## Panel D: interpretability

Program annotations (technical screen, top genes) in `program_annotations.csv` (L0 basis; regenerated on L1-groups when available). Cross-cohort matching against COMBAT (L0 fits, loadings space): **28 matched program pairs across 7 mapped cell types, 13 strong (|cos| ≥ 0.8)** — sex programs replicate in all 6 shared CTs (|cos| 0.979–0.988), interferon programs in B (0.882) and NK (0.880), cytotoxic-CD4 (FCRL6/ADGRG1/GZMH/S1PR5/FGFBP2, 0.846) and monocytic IL1R2/AREG (0.834) programs replicate. See `cross_cohort_matches.csv` and the cross-dataset summary.

## Panel E: practicality

| stage | minutes |
|---|---|
| fit_none (L0, 19 CTs) | 33.9 |
| block_none (L3) | 55.8 |
| fit_groups (L1) | 130.9 |

For reference, the abandoned full-ladder attempt spent >7 h in fit_global without completing — L1-global at 19 CTs is impractical without more wall-clock; the reduced ladder (L0 + L1-groups + block) is the recommended configuration at this scale.

## Limitations

- Reduced ladder: no L1-global, no L2, no PVA/EV-BIDIFAC for this cohort (findings from COMBAT/HLCA carried over: L1-global ≈ L1-groups; L2 unreliable; scores-space integrations degenerate/contaminated).
- Site confounding is strong in this cohort; cp reps reduce but do not eliminate it (cp_L0 0.684 vs floor 0.444).
- Collection_Day/Resample uninformative at donor level.
- Panel D annotations are on the L1-groups fit; block structure is on the L0 basis (block_groups did not finish before the session cap).

## Artifacts

- `eval_all.csv`, `program_summary.csv`, `program_annotations.csv`, `block_structure.csv`, `sharing_spectrum_*.csv`, `runtime.csv` (shared `stephenson/fits/`)
- Figures: `fig1_panelAB_heatmap.png`, `fig2_r2_ladder.png`, `fig3_block_mass.png`
- Fits: `stage1_fits.rds` (L0; groups appended on completion), `integrations.rds` (block_none)
