# cellprograms real-data benchmark: cross-cohort summary

**Scope**: per-cell-type EBMF (flashier) with the local-first borrowing ladder (L0 → L1-global → L1-groups → L2 init; L3 block integration; L4 sharing spectra), benchmarked against 6 cheap sample-level baselines (random, composition, clr, global_pca, perct_pca, grouped_pca) and two comparators (GloScope-GMM, PILOT) on real cohorts. Panels per dataset: A sample-biology preservation, B technical-covariate leakage, C representation quality (R²), D interpretability, E practicality.

| cohort | cells | donors | cell types | ladder | status |
|---|---|---|---|---|---|
| COMBAT (COVID-19 PBMC) | ~1.0M | 50–121/CT | 10 | full (L0/L1g/L1gr/L2 + block/pva/evbS) | **complete** (stability running) |
| HLCA core (lung atlas) | ~584k | 40–69/CT | 9 | full | **complete** |
| Stephenson (COVID-19 PBMC) | ~780k | 120 | 19 | reduced (L0 + L1-groups + block) | **complete** (reduced) |
| OneK1K (PBMC, 300-donor cap) | ~1.3M | 43–300/CT | 9 | reduced (L0 + block) | **complete** (reduced) |

Per-dataset reports: `combat/report_benchmark_combat.md`, `hlca/report_benchmark_hlca.md`, `stephenson/report_benchmark_stephenson.md`, `onek1k/report_benchmark_onek1k.md`.

## Q1 — Do cp representations preserve sample biology?

Consistent pattern across COMBAT and HLCA:

- **cp reps win or tie sex in all four cohorts; win Age in three.** COMBAT: cp_L2 sex 0.859 (best baseline 0.819), cp_L1groups Age 0.464 (baselines ≤0.422). HLCA: cp_L2 sex 0.755 (baselines ≤0.57), cp_L1groups Age 0.322 (baselines ≤0.29; GloScope/PILOT ~0.15). Stephenson: cp_L0 Age 0.653 (baselines ≤0.48, GloScope 0.136) and Days_from_onset 0.559 (baselines ≤0.38); sex a three-way tie at 0.78. OneK1K: **cp_L0 sex 0.978** — the strongest single result in the benchmark. The one Age exception is OneK1K (composition 0.678 wins; cp_L0 0.402 best expression rep) — in a healthy cohort, aging's dominant blood signal is compositional (naive→memory shift).
- **The L3 block representation is the best clinical-label rep on Stephenson**: cp_block wins Outcome (0.861) and ties clr on Status (0.833 vs 0.834), while per-CT programs fragment the systemic COVID signal (cp_L0 Status 0.460). Shared factors capture cross-lineage systemic biology.
- **Compositional phenotypes belong to composition baselines.** COMBAT Death28: composition 0.572 wins. HLCA smoking_status: composition 0.450 / clr 0.392 win. Stephenson Worst_Clinical_Status/Smoker: clr wins (0.475/0.543). When the biology *is* cell-type proportion, expression programs correctly add nothing — a sensible specificity result, not a failure.
- **Cohort/source labels are won by confounded reps.** COMBAT Source: grouped_pca 0.455 (cp ~0.28); HLCA tissue: GloScope 0.754 / perct_pca 0.724 — both study-confounded (Panel B). High "biology" scores that co-vary with batch are not wins.
- **Panel B separation is starkest on Stephenson** (multi-site): perct_pca and grouped_pca hit 1.000 on Site, GloScope 0.960, clr 0.939 (floor 0.444); cp reps are cleanest (cp_block 0.558, cp_L0 0.684).
- **GloScope and PILOT are not competitive as general sample reps here.** GloScope's only wins are confounded; PILOT is near floor on most covariates in all three cohorts.

## Q2 — Are programs stable and interpretable?

- **R² is high and ladder-stable**: COMBAT 0.35–0.85 per CT (rare CTs lowest), HLCA 0.69–0.85.
- **Recurring biology across cell types**: sex programs (XIST/RPS4Y1/DDX3Y/UTY/KDM5D) and interferon programs (IFIT3/IFI44L/RSAD2/USP18/OAS1; cMono instance SIGLEC1-led) recur across COMBAT CTs; HLCA has 2 sex-gene-led + 1 IFN program. Programs are distributed (median top-50 loading mass 0.08–0.10).
- **Per-program technical screen is essential and more sensitive than rep-level kNN.** COMBAT: 2/266 programs Institute η²>0.3 (clean). HLCA: 59/199 study η²>0.3, 37/199 platform η²>0.3 — including the top-PVE program in 5 cell types (η² up to 0.99). Same pipeline, different cohort quality; the screen tells you which programs to distrust.
- **L4 sharing spectra: use loadings space.** Scores-space principal angles degenerate to 0 by construction at ~25–30 programs/CT with ≤~50 shared donors (k_a+k_b > n) on both cohorts. Loadings space is informative: COMBAT cMono↔ncMono 10° (14 shared directions); HLCA AM↔EMac 21.8° (12 shared), RBC most isolated.
- **Cross-cohort replication (COMBAT ↔ Stephenson, independent COVID-19 PBMC cohorts)**: loadings-space matching of L0 fits across 7 mapped cell types recovers **28 program pairs (|cos| ≥ 0.5), 13 strong (≥0.8)**. The strongest matches are biologically exact: sex programs (RPS4Y1/DDX3Y/UTY/EIF1AY/XIST) replicate in all 6 shared CTs (|cos| 0.979–0.988); interferon programs replicate in B (0.882; IFIT3/IFIT1/IFI44L/CMPK2 identical) and NK (0.880); a cytotoxic-CD4 program (FCRL6/ADGRG1/GZMH/S1PR5/FGFBP2) replicates at 0.846; a monocytic IL1R2/AREG program at 0.834. Programs discovered independently per cohort recover the same biology in the same cell types.
- **Stability (COMBAT, 10× subsample refits at 80% donors)**: 132/180 programs (73%) recovered in ≥8/10 refits (median loading cos 0.854); among recovered programs median score cor 0.855. Instability is confined to the low-PVE tail (PVE–recovery Spearman ρ = 0.40; 28 never-recovered programs are minor factors). The interpretable program set is robust to donor resampling.

## Q3 — Is local-first borrowing correct?

- **L1 pooling scope is immaterial at the representation level** in all three cohorts with L1 fits: max rep-level |Δ| between L1-global and L1-groups = 0.039 (COMBAT, sex), 0.027 (HLCA, Age); Stephenson L1-groups ≈ L0 (R² ±0.02 per CT; rep-level within ~0.13 of L0 on all covariates). Program-level matching: COMBAT 251/281 matched (median |cos| 0.967); HLCA 197/199 matched (median |cos| 0.999), identical per-CT program counts. Residual divergence concentrates in **continuous clinical covariates (Age)** of a few programs and is bidirectional. **L1-groups stays the default** (guards cross-lineage homogenization) at no measurable cost — and it is 3.7× (COMBAT) / 2.9× (HLCA) cheaper than L1-global.
- **L2 initialization borrowing is unreliable — negative result replicated.** COMBAT: target rare CTs collapse (PB 25→5 programs, R² 0.57→0.13; GDT 24→8, 0.56→0.20). HLCA: EMac 22→19, R² 0.722→0.651 (RBC mildly helped). Do not use L2 for rare-CT rescue as implemented; L1 rungs are the safe mechanism (HLCA RBC: 11→14 programs, R² 0.810→0.848 under L1).
- **L3 block decomposition is coherent when the cohort is clean.** COMBAT: 9 global + 11 partial factors, lineage-coherent (CD8+DP; cMono+ncMono+DC; CD4+CD8+NK+cMono IFN coalition ×3); classification identical on L0 and L1 bases. HLCA: 16 global + 4 partial, less lineage-coherent — consistent with study-driven sharing. Stephenson (L0 basis): 13 global + 5 partial + 1 private — partials lineage-coherent (two monocyte-subset factors), the single private factor is Platelets (no shared programs with nucleated cells). Block structure is a cohort-quality diagnostic as much as a biology summary.
- **EV-BIDIFAC scores-space states are batch-contaminated in both cohorts** (COMBAT Institute 0.871 vs floor 0.470; HLCA study 0.595/platform 0.508). Its Panel A "wins" are artifacts. Not recommended as a sample representation at this scale.

## Panel E — practicality

| stage | COMBAT (10 CTs, ≤121 donors) | HLCA (9 CTs, ≤69 donors) | Stephenson (19 CTs, 120 donors) | OneK1K (9 CTs, ≤300 donors) |
|---|---|---|---|---|
| fit_none | 20.1 | 15.2 | 33.9 | 16.3 |
| fit_global | 364.4 | 235.5 | >7 h, abandoned | — |
| fit_groups | 98.0 | 81.4 | 130.9 | — |
| fit_l2 | 19.7 | 16.9 | — | — |
| integrations | 75.7 | 76.6 | 55.8 (block only) | 29.7 (block only) |
| **TOTAL** | **578.9** | **426.3** | reduced ladder | **46.4** |

L1-global dominates (~15–18× L0) and becomes impractical at 19 CTs (>7 h); L1-groups delivers equivalent representations for ~3× less. Donor count is not the runtime driver (OneK1K at 300 donors is the fastest) — L1 prior iterations are. GloScope-GMM runs in minutes (KNN densitometer produced pervasive non-finite distances on COMBAT — dropped); PILOT is fast but weak. Block integration ~30–77 min.

## Recommendations (current evidence)

1. Default pipeline: **L1-groups + block integration**, loadings-space L4; skip L1-global unless testing pooling scope; **do not use L2** for rare-CT rescue as implemented.
2. Always run the **per-program technical screen**; rep-level Panel B alone misses program-level contamination (HLCA).
3. Interpret "biology wins" by composition/GloScope only after checking the matching technical covariate.
4. At ≥25 programs/CT and ≤~120 donors, report L4 in loadings space only.

## Caveats

- Four cohorts evaluated; Stephenson and OneK1K ran reduced ladders (L0 + L1-groups + block / L0 + block) under the session time cap — their L1-global/L2/PVA/EV-BIDIFAC cells are covered by the COMBAT/HLCA findings.
- OneK1K has low reconstruction R² (0.16–0.36): healthy cohort, 300 donors, capped factor count. Its programs are nonetheless biologically real (sex 0.978).
- OneK1K pool_number per-program screen is uninformative (73 pools / 300 donors inflates η²; null ≈ 0.24).
- kNN-in-266-dim sparse program space vs 10–30-dim PCA is a metric×dimensionality comparison; cp wins on sex/Age are conservative under this handicap, but absolute values are not directly comparable across representation families.
- HLCA is multi-study; its block/global structure partly reflects batch.
