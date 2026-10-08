# HLCA core benchmark: cellprograms vs sample-level baselines and comparators

**Dataset**: Human Lung Cell Atlas (core, ~584k cells) → donor × cell-type pseudobulk, 9 cell types (AM, EMac, cMono, DC, CD8, MCC, AT2, RBC, VEC), 40–69 donors per cell type, top-3000 HVGs. Multi-study atlas — *study* and *sequencing_platform* are the technical covariates; biology covariates are smoking_status, sex, tissue_level_2 (categorical) and Age (continuous).
**Method**: per-cell-type EBMF (flashier), borrowing ladder L0 → L1-global → L1-groups (lineage groups: Myeloid/Epithelial/T/Endothelial) → L2 init; L3 block integration (20 shared factors, expression level); L4 sharing spectra. Comparators: 6 cheap baselines + GloScope (GMM) + PILOT. Config: max_factors=25, prior_max_iter=3, seed 42.

## Headline results

1. **cp representations carry the strongest sex and Age signal.** cp_L2 wins sex (kNN macro-F1 0.755 vs ≤0.57 for all baselines and comparators); cp_L1groups wins Age (max|Spearman| 0.322 vs ≤0.29 baselines; GloScope 0.157, PILOT 0.147). Same pattern as COMBAT.
2. **Composition/clr win smoking_status (0.39–0.45)** — a genuinely compositional phenotype; expression programs add little. Consistent with COMBAT's Death28 result.
3. **GloScope "wins" tissue_level_2 (0.754) but is study-confounded (0.810)** — tissue and study are entangled in HLCA; perct_pca shows the same (0.724/0.756). cp reps keep both low (tissue ≤0.28, study ≤0.41).
4. **Panel B: cp reps are the cleanest expression representations** (study 0.354–0.406, platform 0.286; random floor 0.093/0.069) vs global_pca 0.722/0.578, grouped_pca 0.623/0.602, GloScope 0.810/0.459. Only cp_block (0.314/0.126) and PILOT (0.062/0.118, but near-floor on biology too) are lower.
5. **The per-program technical screen is more sensitive than rep-level kNN**: 59/199 programs have study η²>0.3 and 37/199 platform η²>0.3 — including the *top-PVE* program in several cell types (CD8__F001 η²=0.99, VEC__F001 0.99, AT2__F001 0.98). HLCA programs must be interpreted study-aware; the screen flags exactly which ones to distrust.

## Panel A/B: covariate preservation

![fig1](fig1_panelAB_heatmap.png)

kNN macro-F1 (categorical) / max |Spearman| (continuous). Panel A left of the divider, Panel B right.

- **sex**: cp_L2 0.755 > cp_L0/L1 0.657–0.674 > global_pca 0.567 > GloScope 0.541 > PILOT 0.478. Only 2/199 programs are explicitly sex-gene-led (MCC__F023, AT2__F021; XIST/RPS4Y1/DDX3Y), both low-PVE — the sex signal is distributed across many programs, and the representation still wins.
- **Age**: cp_L1groups 0.322, all cp variants ≥0.295; baselines 0.26–0.29; GloScope/PILOT ~0.15.
- **smoking_status**: composition 0.450 / clr 0.392 win; cp_L1global 0.360 best expression rep. Smoking's signature here is cell-compositional.
- **tissue_level_2**: GloScope 0.754 and perct_pca 0.724 "win" — but both are study-confounded (see Panel B); tissue is not cleanly separable from study in this atlas.
- **Panel B**: cp_evbS is again contaminated (study 0.595, platform 0.508, tissue 0.746) — second cohort where EV-BIDIFAC scores-space states absorb technical structure. cp_block is the cleanest cp variant (platform 0.126).

## L1-global vs L1-groups

Rep-level associations are nearly identical: max |Δ| per covariate = 0.027 (Age), 0.018 (sex), ≤0.004 elsewhere. Program-level matching confirms it: identical per-CT program counts, **197/199 programs matched (median |cos| 0.999)**, 4 unmatched rows are ±1 cap swaps in AM/cMono. As in COMBAT, the residual divergence concentrates in **Age associations** of a few programs (max |Δ| 0.345, VEC__F015: 0.463 global vs 0.118 groups) and is bidirectional (AM__F005: 0.006 vs 0.218). Second cohort in which the L1 pooling scope is immaterial at the representation level; L1-groups remains the default on principle (guards against cross-lineage homogenization) at no measurable cost.

## Panel C: representation quality and the L3 view

![fig2](fig2_r2_ladder.png)

Per-cell-type R² is high and stable across the ladder (0.69–0.85). Borrowing effects concentrate in the lowest-coverage cell type: **RBC (40 donors) gains 11→14 programs under L1 (R² 0.810→0.848)** — the intended rare-cell-type rescue. **L2 again hurts a target**: EMac 22→19 programs, R² 0.722→0.651 (RBC is mildly helped: 11→12, 0.810→0.818). Net: L2 initialization borrowing is unreliable — consistent with the COMBAT collapse (PB 25→5, GDT 24→8 programs).

![fig3](fig3_block_mass.png)

Block decomposition (L1-global basis): **16/20 shared factors are global, 4 partial** (MCC+VEC, RBC-only, EMac+RBC, AM+EMac+CD8+AT2). Less lineage-coherent than COMBAT — plausible given the study-driven program structure above: part of what is "shared" across HLCA cell types is study batch structure, not lineage biology.

## Panel D: interpretability

- 199 programs; median top-50 loading mass 0.081 — distributed programs, as in COMBAT (0.10).
- Loadings-space sharing spectrum is informative: AM↔EMac are the closest pair (min angle 21.8°, 12 shared directions) — biologically coherent (two macrophage compartments); RBC pairs are the most isolated (34–39°). Scores-space spectra are degenerate at this program density (k_a+k_b > n_shared donors), as on COMBAT — use loadings space.
- Technical screen: 59/199 programs study η²>0.3; 37/199 platform η²>0.3. The highest-PVE program in CD8, VEC, AT2, AM, EMac is near-perfectly study-confounded — **do not interpret top programs in HLCA without consulting the screen**.

## Panel E: practicality

| stage | minutes |
|---|---|
| fit_none (L0) | 15.2 |
| fit_global (L1) | 235.5 |
| fit_groups (L1) | 81.4 |
| fit_l2 | 16.9 |
| integrations (block×2 + pva + evbS) | 76.6 |
| **TOTAL** | **426.3** |

L1-global is again the dominant cost (~15× L0). GloScope-GMM ran in minutes; PILOT is fast but near-floor on every HLCA covariate (≤0.22 except sex 0.478).

## Limitations

- HLCA is multi-study: study/tissue confounding limits clean interpretation of Panel A "wins" by GloScope/perct_pca and of the block decomposition's global factors.
- Scores-space L4 spectra degenerate (0/36 pairs significant) — loadings space recommended at ~25 programs/CT and ≤69 shared donors.
- No stability (subsample) analysis for HLCA in this run; COMBAT stability is the template.

## Artifacts

- `eval_all.csv`, `program_summary.csv`, `program_annotations.csv`, `block_structure.csv`, `sharing_spectrum_*.csv`, `runtime.csv` (shared `hlca/fits/`)
- Figures: `fig1_panelAB_heatmap.png`, `fig2_r2_ladder.png`, `fig3_block_mass.png`
- Fits: `stage1_fits.rds`, `integrations.rds` (shared `hlca/fits/`)
