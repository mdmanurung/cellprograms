# Coordinated multicellular program (CMP) detection in `cellprograms`: v2 report

Project: CellPrograms. Repository: `mdmanurung/cellprograms`. Workflow: tusoskill benchmark-first iteration (58 logged substantive iterations; floor 50).

## 1. Summary

- **What changed.** The CMP caller was already frozen after v1.1 (`n_agg_noanc`: rank-first scores, ambient filter, node-aggregate max-T). Plan v2 tested three routes to better detection: better upstream programs, joint factorization (MOFA+), and caller changes (breadth-adaptive and context-supervised tests). Only the upstream route helped: **refitting per-cell-type programs with `max_factors = 30` instead of 10**.
- **Held-out confirmation (pre-registered cmp-bench-v1.2, seeds 1051-1100, never used before).** All 7 criteria passed.
  - Pooled membership F1 0.705 vs 0.571 for the v1.1 finalist (paired difference +0.127, 95% CI 0.113-0.142, n = 600 positive datasets).
  - The gain is entirely on the semi-synthetic base built from Stephenson pseudobulk (F1 0.671 -> 0.934, +0.251). The parametric base is neutral (+0.003), because its simulated cell types rarely carry more than 10 programs (maximum observed 14).
  - Versus repo baselines: +0.353 over node correlation and +0.351 over principal angles.
  - False-call guardrails (parametric / semi): null FWER 0.04 / 0.00; ambient decoy 0.04 / 0.02; within-context false 0 / 0; genes-only false 0.06 / 0. Shuffled control F1 0.0003.
- **Package.** The global `max_factors` default stays at 10 (deviation from the v1.2 pre-registration wording; see section 6). `find_coordinated_programs()` documents the cap-30 recommendation and prints a one-time message when a fit with a lower cap has a cell type at its cap. Tests added for LR FDR control, context-association missingness, and the message. Full suite: 42 tests, 0 failures.
- **Real data.** Refitting Stephenson and COMBAT at cap 30 roughly doubles the CMP count (13 -> 27; 11 -> 25). The CMPs kept from cap 10 retain their annotations, and cross-cohort matching improves (4 -> 7 matched CMPs). Many of the 13 new CMPs per cohort carry immunoglobulin/TCR V genes or genes with known donor-genotype variation (GSTM1 deletion, HLA-DQ/DR haplotype, CCL3L3/CCL4L2 copy number), and only 1 of 26 is associated with severity (p < 0.05). How many count as V(D)J- or genotype-marked depends on an unvalidated heuristic: 7/13 (Stephenson) and 9/13 (COMBAT) with >= 1 marker gene, but 3/13 and 8/13 with >= 2. A genotype origin is therefore a hypothesis to test with genotype or eQTL data, not an established finding. Cap 30 is still preferable to 40 on real data, and real-data CMPs from high-cap fits need annotation-based triage (section 5).

## 2. Benchmark design (unchanged from v1.1)

- **Simulation.** 8 scenarios: null, cmp_distinct (CMP members with distinct genes), cmp_mixed, context (CMP driven by a donor context variable), ambient (ambient-RNA decoy), lr (ligand-receptor-linked CMP), weak (low amplitude), genes_only (shared genes without donor co-variation, a decoy). Each positive dataset plants CMPs of different breadth: 20%, 50% and 100% of 12 cell types, plus a context-restricted CMP.
- **Two bases.** Parametric (simulated programs, amplitude 0.45) and semi-synthetic (real Stephenson COVID-19 PBMC pseudobulk, 12 major cell types, with planted donor-level signal; amplitude 0.875).
- **Primary endpoint.** Membership F1 of called CMPs against planted members, averaged over positive scenarios (marginal mode).
- **Upstream recall ("oracle").** Best F1 achievable given the fitted programs, i.e. whether a program matching each planted member exists at all.
- **Guardrails.** Null FWER, ambient decoy call rate, within-context false calls, genes-only false calls, each <= 0.10 per base.
- **Splits.** Dev seeds 1-30; null calibration seeds 31-130; v1.1 test 1001-1050; v1.2 held-out 1051-1100.

## 3. What was tried (plan v2)

| Route | Arm | Result on dev | Decision |
|---|---|---|---|
| Upstream | `max_factors` 20 (seeds 1-15) | +0.115 [0.088, 0.143]; semi +0.228, param +0.001 | promoted to cap screen |
| Upstream | `max_factors` 30 / 40 (semi seeds 1-8) | F1 0.9325 for both vs 0.690 at K10; upstream recall 1.00 | K30 chosen (40 is 23% slower, no gain) |
| Upstream | `max_factors` 30, full dev 1-30 + nullx | 0.692 vs 0.577; +0.111 [0.092, 0.130] | to held-out |
| Upstream | all genes / 4000 HVGs (param) | -0.211 / -0.218; ambient decoy 0.40 / 0.33 | rejected |
| Upstream | point-normal prior | -0.091 [-0.120, -0.064] | rejected |
| Joint | MOFA+ cross-fitted hybrid (J-A) | F1 0.047; param null FWER 0.6 (n = 5) | rejected |
| Joint | MOFA+ augment (J-B) | -0.007 [-0.023, 0.006] | rejected |
| Caller | C1 breadth-adaptive | +0.002 [-0.002, 0.005]; semi genes-only false 0.10 | not adopted |
| Caller | C2 context-supervised | -0.002 [-0.011, 0.007] | not adopted |
| Caller | C1 + C2 | +0.004 [-0.002, 0.011] | not adopted |

Variants that reached or broke a guardrail, with their F1 change, are collected in `tusoskill/cmp_v2/tradeoffs.csv`. None was adopted. The two v1.1 runners-up (`n_agg_t2`, edge max-T) gained at most +0.004 F1 on the v1.1 test set while raising parametric null FWER to 0.14.

## 4. Held-out results (cmp-bench-v1.2)

Figures are in `tusoskill/cmp_v2/figures/`:

- `fig1_f1_by_scenario_heldout.png`: F1 by scenario and arm.
- `fig2_miss_by_breadth_heldout.png`: miss rate by planted CMP breadth.
- `fig3_upstream_recall_heldout.png`: upstream recall, cap 10 vs cap 30.
- `fig4_null_fwer_heldout.png`: null FWER with Clopper-Pearson CIs (held-out nulls, and pooled with 100 extra null datasets per base).
- `fig5_paired_forest_heldout.png`: paired differences, winner minus each comparator.

| Run | Pooled F1 | Parametric F1 | Semi F1 |
|---|---|---|---|
| Cap 30 + `n_agg_noanc` (winner) | 0.705 | 0.475 | 0.934 |
| Cap 10 + `n_agg_noanc` (v1.1 finalist) | 0.571 | 0.472 | 0.671 |
| Node correlation baseline | 0.349 | 0.018 | 0.679 |
| Principal angles baseline | 0.355 | 0.310 | 0.400 |
| Shuffled control (cap 30) | 0.0003 | | |

| Comparator | Paired diff [95% CI] | Parametric | Semi |
|---|---|---|---|
| v1.1 finalist | +0.127 [0.113, 0.142] | +0.003 [0.001, 0.006] | +0.251 [0.229, 0.273] |
| Node correlation | +0.353 [0.332, 0.375] | +0.463 | +0.243 |
| Principal angles | +0.351 [0.328, 0.373] | +0.165 | +0.537 |

**Mechanism.** On the semi base, cap 30 raises upstream recall from 0.73-0.83 to 0.995-1.00 across breadths (fig 3). Narrow-CMP (20% of cell types) misses fall from 0.55 to 0.16, and context-restricted misses from 0.08 to 0. Real cell types carry more than 10 identifiable programs; at cap 10, the CMP member is often never fitted.

**What did not improve.** On the parametric base, narrow CMPs are still missed 71% of the time and the context-restricted CMP 66% of the time, at both caps. The principal-angles baseline misses fewer narrow parametric CMPs (61%), but it calls an ambient decoy in 100% of ambient datasets. Detecting a 2-of-12 cell-type CMP at amplitude 0.45 appears power-limited, and no caller variant changed that.

**Null calibration.** Pooling held-out and extra nulls (150 per base), the winner's FWER is 0.033 (parametric; CP 95% CI 0.011-0.076) and 0.000 (semi; 0-0.024).

## 5. Real data: Stephenson and COMBAT refitted at cap 30

**Setup.** Programs refitted with `max_factors = 30` and covariates institute + sex. The frozen caller was run in flag mode with 999 permutations. Compared with the cap-10 runs (Stephenson with sex covariate; COMBAT without) and cap-40 runs (both with sex). Cap-to-cap matching uses |Pearson r| between CMP consensus donor scores on the same donors; a CMP is "retained" if its best match is >= 0.5, otherwise "new" (or "lost" for the reference run). The table is `tusoskill/cmp_v2/real_K10_K30_K40_summary.csv`, with per-CMP tables in `tusoskill/cmp_v1/real/*_vs_*`.

Annotation flags (computed on each CMP's 12 most recurrent top genes):
- V(D)J-dominated: >= 50% match `^(IG[HKL]V|IG[KL]C|IGH[GAMDE]|JCHAIN|TR[ABGD]V)`.
- Genotype-marked: >= 1 gene from a fixed list of genes with common deletion, copy-number, HLA-haplotype or sex-chromosome expression variation (GSTM1, GSTT1/2, UTS2, ERAP2, CCL3L1/3L3/4L2, HLA-DRB3/4/5, HLA-DQA1/2, HLA-DQB1/2, XIST and Y genes). This is a screening heuristic, not validated against genotype or eQTL data. A single marker gene can tag a CMP that is otherwise dominated by other biology (e.g. COMBAT CMP23, complement/C1Q genes plus CCL3L3). With a stricter >= 2 marker genes, the genotype-marked counts at cap 10 / 30 / 40 are 0 / 4 / 4 (Stephenson) and 1 / 4 / 2 (COMBAT).

| Cohort | Cap | CMPs | >= 5 cell types | 2-4 cell types | V(D)J-dominated | Mostly ambient | Severity p < 1e-3 | Genotype-marked | New vs cap 10 | Lost vs cap 10 | Cross-cohort matches |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Stephenson | 10 | 13 | 11 | 2 | 4 | 4 | 3 | 1 | | | 4 |
| Stephenson | 30 | 27 | 22 | 5 | 2 | 4 | 2 | 8 | 13 | 2 | 7 |
| Stephenson | 40 | 31 | 24 | 7 | 6 | 3 | 5 | 9 | 18 | 3 | 5 |
| COMBAT | 10 | 11 | 8 | 3 | 1 | 4 | 3 | 2 | | | 4 |
| COMBAT | 30 | 25 | 15 | 10 | 6 | 4 | 5 | 6 | 13 | 2 | 7 |
| COMBAT | 40 | 30 | 15 | 15 | 9 | 4 | 3 | 11 | 18 | 1 | 5 |

Cross-cohort matches count Stephenson CMPs whose best COMBAT match has mean |loading r| >= 0.5 over >= 3 shared cell types. The COMBAT cap-10 run lacks the sex covariate, so its cap-10 vs cap-30 comparison changes two things at once.

**Retained CMPs.** All 11 multi-cell-type CMPs from Stephenson cap 10 are retained at cap 30 (|r| 0.67-0.99):
- the pan-cell-type interferon CMP (severity p 0.042 -> 0.008);
- the ambient immunoglobulin CMP (p = 7e-5);
- proliferation (p = 3e-15);
- the platelet, monocyte-ambient, cytotoxic and TNF/NF-kB CMPs.

One cap-10 severity hit weakened: CMP10 (allograft/IL2-STAT5, p = 2e-4) maps to cap-30 CMP21 at p = 0.065. In COMBAT, 9 of 11 are retained. The two lost CMPs are an XIST-containing sex signal (removed by the sex covariate, as expected) and an interferon-gamma/inflammatory CMP that maps partly (|r| = 0.50).

**Lost CMPs.** Stephenson loses only the two 2-cell-type V(D)J CMPs at cap 30.

**New CMPs (13 per cohort).**
- Stephenson: 6 genotype-marked (2 with >= 2 marker genes); 4 narrow (2-4 cell types); 1 V(D)J-dominated by the >= 50% rule, though most others also carry several V genes.
- COMBAT: 6 genotype-marked (4 with >= 2 marker genes); 7 narrow; 4 V(D)J-dominated; 1 an immediate-early stress signature (JUN/FOS/ATF3).
- Several other new CMPs carry genes with known strong cis-eQTL or copy-number variation that are not in the fixed list (RPS26, GATD3, POMZP3, CCDC144A).
- One new Stephenson CMP (CMP25, 2 cell types) contains RPS4Y1, USP9Y and DDX3Y despite the sex covariate. This is worth checking against the donor sex labels.
- Only 1 of the 26 new CMPs is associated with severity at p < 0.05: Stephenson CMP19, an 11-cell-type interferon-flavoured CMP with 4 of 11 members ambient-flagged (p = 0.003).

**Interpretation.** Donor genotype produces real cross-cell-type co-variation: a donor homozygous for the GSTM1 deletion lacks GSTM1 in every cell type. If that is the source, the caller is right to call such CMPs, but they are not coordinated biology. This is a hypothesis: no genotype data were used, and the marker heuristic is sensitive to its threshold (above). Consistent with a genotype origin, the GSTM1 CMP is one of the cross-cohort matches at cap 30 (Stephenson S9 vs COMBAT C13, r = 0.54 over 6 cell types). The other robust matches are interferon, immunoglobulin, monocyte (S100A8/9), AREG/SOCS3 cytokine response, proliferation and an IFI27-containing CMP.

**Site replication (Stephenson; Ncl 62 and Cambridge 46 donors, up to 60 / 45 shared per cell-type pair; programs learned on all donors).**

| Cap | Edge sign concordant, Cambridge / Ncl | Concordant and nominal p < 0.05, Cambridge / Ncl |
|---|---|---|
| 10 | 0.922 / 0.977 | 0.621 / 0.770 |
| 30 | 0.897 / 0.958 | 0.548 / 0.740 |
| 40 | 0.902 / 0.971 | 0.542 / 0.771 |

The new narrow Stephenson CMPs (17, 20, 25, 27; 2-3 cell types each) have edge signs concordant in both sites (1.0 / 1.0 for three of them; 0.67 / 1.0 for CMP27). Nominal replication in Cambridge, the smaller site, varies: 0, 0.67, 1.0 and 0.67. Independent re-detection within a site recovers at most 2 of their members (CMP17: 2 in Ncl, 0 in Cambridge; CMP20: 2 / 0; CMP25: 0 / 0; CMP27: 0 / 1). They are reproducible in direction but underpowered within a single site, and none is associated with severity. Treat them as candidates.

**Cap 30 vs 40 on real data.** Every COMBAT cap-30 CMP has a cap-40 match (0 new). The 5 cap-40 CMPs absent at cap 30 are 2-3-cell-type groups with no severity association. In Stephenson, cap 40 has 8 CMPs without a cap-30 match. One is a 6-cell-type interferon CMP (severity p = 6e-4) that partly maps to cap-30 CMPs (best |r| = 0.48); the rest are narrow or have no severity association. Cap 40 adds more V(D)J-dominated CMPs (6 vs 2 in Stephenson, 9 vs 6 in COMBAT) and fewer cross-cohort matches (5 vs 7). This supports 30 over 40, consistent with the simulation finding that 40 adds fit time without recall gain.

**Practical recommendation.** For real cohorts, fit at cap 30 and call in flag mode. Then triage each CMP using:
- ambient flags;
- V(D)J dominance;
- genotype markers (and, where available, genotype data or known eQTL lists);
- severity/context association;
- site replication.

This step is needed because the simulation benchmark contains no V(D)J or genotype structure (section 7).

### 5.1 TF footprints and triage of the cap-30 CMPs

Applied `annotate_cmp_tf()` (CollecTRI from OmniPath: 42,601 signed TF-target pairs after merging NF-kB/AP-1 complexes; TFs with >= 5 targets among each program's 2,000 genes; median 300 TFs tested per member) and `cmp_triage()` (severity association adjusted for institute and sex; Stephenson sign replication within Cambridge and Ncl; COMBAT has no second site with >= 15 donors). Outputs: `tusoskill/cmp_v2/interpret/`.

| | Stephenson (27 CMPs) | COMBAT (25 CMPs) |
|---|---|---|
| candidate / candidate, severity p < 0.05 | 10 / 6 | 6 / 3 |
| ambient | 4 | 6 |
| V(D)J | 1 | 5 |
| genotype candidate (>= 2 marker genes) | 3 | 4 |
| platelet / erythrocyte | 2 | 0 |
| sex chromosome | 1 | 0 |
| stress (IEG) | 0 | 1 |
| low within-site replication (soft note) | 7 | n/a |

TF footprints:
- **IFN CMPs show STAT1/STAT2/IRF9 footprints in both cohorts** (Stephenson CMP1: STAT1 in 10/11 members; COMBAT CMP2: 9/9). This is not an independent positive control: the competitive p-value was switched from a gene-permutation p (floored at 1/1001, which capped BH and called STAT1 in only 5/11 Stephenson CMP1 members, with no shared TF) to the parametric ULM t-test p after these real-data results had been seen. The parametric p assumes independent genes and is more anticonservative for co-regulated blocks.
- **The ambient monocyte CMP (S3) carries a CEBP/SPI1 footprint in 8/10 members, 7 of its 9 non-myeloid members** (the eighth is DC, which is myeloid), consistent with monocyte transcript leakage into T, B and NK programs.
- **COMBAT CMP1 (severity p = 3e-4) has a shared RELA/FOS/ATF2/CREB1 footprint** (NF-kB/AP-1 activation with AREG, NFKBIA, CD83).
- **Stephenson CMP7, labelled genotype candidate (GSTM1 plus V genes), has a strong negative CREB1/ELK1/EGR1/RELA/ATF2 footprint in 8/8 members**, i.e. an immediate-early axis; its footprint anti-correlates with COMBAT CMP17 (stress_ieg) at r = -0.81. The genotype label is likely wrong for this CMP: the marker genes ride along with a stress axis.
- **COMBAT CMP24 (HLA-DQA1/DRB5, genotype candidate) shows RFX5/RFXAP/CIITA.** This is circular, because the HLA class II genes are those TFs' targets, so it neither supports nor refutes a genotype origin.
- **Cross-cohort consistency.** For the 7 robust Stephenson-COMBAT matches (section 5), the matched partner has the highest |footprint correlation| among all 25 COMBAT CMPs in every case (|r| 0.59-0.95; S19-C18 is anti-aligned, r = -0.63), vs a 95th percentile of 0.41 over all pairs. Matched CMPs share genes, so this shows reproducibility of footprints, not independent validation.

Caveats: the competitive test treats genes as exchangeable, so co-regulated gene blocks inflate significance; large regulons (NFKB, AP1, SP1) are called more readily; footprints are descriptive, not causal. The triage labels are heuristic, and the genotype list is unvalidated (section 7).

## 6. Deviations from the pre-registration and plan

1. **Package default kept at 10.** The v1.2 pre-registration says all criteria "must hold to adopt cap 30 as the package default". All held, but the global default is not changed. An earlier repository benchmark (32-config factorial on COMBAT pseudobulk, disease-status kNN and severity endpoints) found that `max_factors >= 30` overfits the concatenated representation (400-600 factors) and collapses disease-signal recovery. The two benchmarks optimize different targets, so the recommendation is scoped to CMP detection: documentation plus a one-time message.
2. **v1.1 test set reduced to 50 datasets per base** (recorded at the time).
3. **Write collisions.** Three weak-scenario held-out fits (seeds 1093-1095) were written by two machines. The existing files were verified intact; fits are deterministic (seed 1).
4. **Sandbox expiry.** 112 semi cap-30 held-out fits were lost in a 24-hour sandbox expiry and refitted with identical code and seeds before scoring. Each was checked for seed, scenario, arm and amplitude.
5. **Bootstrap resamples.** The v1.2 checker uses 2000 paired resamples, as pre-registered (the v1.1 checker used 5000).

## 7. Limitations and what I would change

- **No V(D)J or genotype simulation.** At cap 30, many new real-data CMPs carry immunoglobulin/TCR V genes or candidate donor-genotype marker genes (section 5; 7/13 and 9/13 at >= 1 marker gene, 3/13 and 8/13 at >= 2). Whether they arise from genotype is untested. The simulator has neither clonotype-like nor genotype-like donor-specific gene blocks, so their rate in calls is unmeasured, and the guardrails do not cover them. Next steps: V(D)J and germline-variant decoy scenarios, and an optional genotype-aware filter (e.g. regress out known cis-eQTL genotypes, or flag CMPs whose top genes are dominated by eQTL/CNV genes).
- **Narrow amplitude basis.** Amplitudes were tuned once per base (0.45 / 0.875). Power curves over amplitude and donor number would show where narrow CMPs become detectable.
- **The cap issue was found late.** Upstream recall should have been a first-class diagnostic from v1.0; the K = 10 cap explained most of the semi-base gap.
- **No simulated replication scenario.** Site replication is assessed only on real data; a simulated two-site design would calibrate what concordance to expect.
- **Bases build ambient profiles differently.** Parametric ambient profiles are simulated; semi-base profiles come from real pseudobulk. Ambient guardrails are therefore not strictly comparable across bases.
- **Real data has no ground truth.** New CMPs at cap 30 are candidates; replication and annotation evidence is reported, not proof.

## 8. Files

- `tusoskill/cmp_v2/heldout_v1.2_scorecard.csv`, `heldout_v1.2_paired.csv`, `heldout_v1.2_criteria.csv`: held-out scoring.
- `tusoskill/cmp_v2/tradeoffs.csv`: guardrail trade-offs.
- `tusoskill/cmp_v2/screen_*.csv`: v2 screens.
- `tusoskill/cmp_v2/figures/`: figures and their source tables.
- `tusoskill/cmp_v2/bench_scripts/`: all benchmark scripts and pre-registrations.
- `tusoskill/cmp_v1/real/`: real-data outputs per run, cap comparisons and cross-cohort matching.
