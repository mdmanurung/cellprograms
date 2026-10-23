# `cellprograms` gene-sharing backend: revised implementation and benchmarking plan (v2)

**Date:** 2026-10-10
**Supersedes:** "Audited Hierarchical EBMF Implementation and Benchmarking Plan" (2026-10-09)
**Status:** Proposed. A toy feasibility test (Section 6) has been run. Nothing has been implemented in the package. Repository claims are not yet re-verified (Section 3).

---

## 0. What changed from v1, and why

| # | v1 position | v2 position | Reason |
|---|---|---|---|
| 1 | Build a custom coupled optimizer (M3) and add empirical Bayes only later (M6) | **Fit the exact-core model as `flashier` on row-stacked views, with a custom block-wise score prior.** Use the same EBMF objective (ELBO) from the start | Mathematically equivalent model (Section 5.1). The toy test shows it runs and borrows strength (Section 6) |
| 2 | Membership from the fitted model (h_ck) | **Membership from a cross-fitted, multiplicity-corrected per-view test.** Fitted variances v_b > 0 are *not* membership calls | Toy test: naive calls gave 65–75% false sharing on null data. The cross-fitted test gave 0/40 |
| 3 | Competitors: B0, B1, LIGER, MOFAcellulaR | **Add MOFA+ multi-group, B1 + residual recruitment, and MSFA/Tetris (small grid).** MOFAcellulaR moves to post-hoc coordination only | The closest prior art was missing. MOFAcellulaR estimates a different quantity (coordinated donor scores) |
| 4 | Signal strength as raw SNR | **Signal strength relative to the single-view and pooled detection thresholds** | Borrowing helps in a predictable range. This makes the design efficient and interpretable |
| 5 | "Meaningful improvement" left undefined | **One preregistered primary endpoint with a margin.** Methods compared at a matched null false-sharing rate | Removes freedom to choose the analysis after seeing results |
| 6 | Shared loading vector defined loosely | **Preprocessing contract written first** (centring, scaling, gene-eligibility mask, precision weights) | Without it, "exact signed core" is not well defined |
| 7 | Tier 3 semi-synthetic data used as the realistic null | **Add a stratified within-view gene-permutation null** and real-data positive/negative controls | Gives a null that matches each view's own spectrum and noise structure |
| 8 | Full G0–G9 × 3 tiers × 8 methods grid | **A compact pilot comes first; the preregistration is frozen after the pilot** | The large grid only makes sense once the pilot shows joint fitting adds something |

---

## 1. Objective, estimands and definitions

### 1.1 Objective
Discover **fully shared, partially shared and private** transcriptional programs across cell-type views of pseudobulk scRNA-seq data. The method must not require matched donors or coordinated sample activity. It must also allow the outcome "no shared programs". Donor coordination and phenotype association are analysed **post hoc** only.

### 1.2 Data
For view (cell type) c = 1..C there is a pseudobulk matrix X_c (n_c × G). Rows are pseudobulk observations (donor, or donor × sample). The views share a gene universe 𝒢. Observation sets may differ, overlap or be disjoint across views.

### 1.3 Primary estimand (v1)
A program k is a unit-norm signed gene-loading vector g_k ∈ ℝ^G, defined up to a global sign. Its **membership set** A_k ⊆ {1..C} contains the views where the term z_ck g_kᵀ contributes variance above zero. Program classes:

- **Shared:** |A_k| ≥ 2 (fully shared if A_k equals all eligible views).
- **Private:** |A_k| = 1.
- **Indeterminate:** the test for a view has too little power to call it (for example, too few observations or genes eligible in that view). Report this; never count it as absent.
- **Unobserved:** the view lacks the program's genes or was not fitted.

### 1.4 Secondary quantities (reported, not used for selection)
| Quantity | Measures |
|---|---|
| Gene-support overlap (Jaccard index / AUPRC) | Which genes contribute |
| Loading-subspace similarity (principal angles) | Robustness to factors splitting or rotating |
| Per-view loading deviation ‖w_ck/‖w_ck‖ − g_k‖ | Approximate conservation (M5 onwards) |
| Sample-score coordination | Post hoc, matched donors only |

### 1.5 Scientific hypotheses (preregistered)
- **H1 (borrowing):** In the "borrowing window" (Section 7.2), joint fitting recovers shared programs with correct membership more often than independent EBMF + consensus + residual recruitment (B1+R), at a matched null false-sharing rate.
- **H2 (no harm):** Joint fitting is non-inferior to independent EBMF (B0) for private-program recovery.
- **H3 (calibration):** On structural nulls, the family-wise false-sharing rate is ≤ 5% (upper 95% CI ≤ 8%).

H1 is the main claim. If H1 fails, development stops at M4 (Section 11).

---

## 2. Prior art and positioning

| Method | Model | Relation to this plan |
|---|---|---|
| **MOFA+ multi-group** (Argelaguet et al. 2020) | Weights shared across groups. Each group has its own factor values. Group-level sparsity lets a factor be inactive in a group. Features are centred per group, and mean differences between groups are not modelled | **Closest published competitor** when cell types are the *groups* in one view. Does not need matched samples |
| **MSFA / BMSFA** (De Vito et al. 2021) | Factors common to all studies plus study-specific factors | Fully shared vs private, but no partial sharing |
| **Tetris** (Grabski et al. 2023) | An Indian Buffet Process over the subsets of studies that share each factor | **The same estimand as this plan**, solved by Gibbs sampling. Practical only on small grids |
| **SLIDE** (Gaynanova & Li 2019), **JIVE** (Lock et al. 2013) | Joint/individual or partially shared structure | Similar ideas, but sharing runs across samples (views share samples, not features). The axis is transposed |
| **LIGER iNMF / UINMF** (Welch et al. 2019; Kriebel & Welch 2022) | Shared W plus dataset-specific V_i, non-negative. UINMF allows features unique to one dataset | The structure behind the v2 flexible model (M5). Non-negativity is the main difference |
| **Spectra** (Kunes et al. 2023) | Global and cell-type-specific gene programs, guided by gene sets | Supervised. Useful as a real-data reference, not a primary competitor |
| **cNMF, GeneNMF** | Consensus programs across runs or samples | Basis of the post-hoc consensus baseline B1 |
| **MOFAcellulaR** (Ramirez Flores et al. 2023), **DIALOGUE** (Jerby-Arnon & Regev 2022) | Coordinated multicellular programs (shared donor scores) | **Different estimand.** Used only in the post-hoc coordination analysis |

**Novelty claim, if H1–H3 hold:** an empirical Bayes factorization with sparse, adaptive gene-loading priors that:
- infers which subsets of cell-type views share each program, with calibrated membership tests;
- runs at atlas scale without matched donors;
- needs no MCMC (unlike Tetris);
- has no non-negativity constraint (unlike LIGER);
- and reduces exactly to independent EBMF when nothing is shared.

Claiming "the first method for partially shared programs" would be wrong, because Tetris already does this. The claim should be about **calibrated membership, scale, and adaptive sparse loadings**.

---

## 3. Repository audit: verification checklist

The repository was not readable from this session (GitHub API returned 404, probably because it is private). So the v1 audit table is **unverified**. Each item becomes a test that must fail on current `main` before its fix is merged:

| ID | Claim (v1) | Verification test (write first) | Fix |
|---|---|---|---|
| A1 | `mofacellular.R` combines matrices with `cbind` without aligning genes | Two views with reversed gene order and one view with an extra gene. Expect alignment by name or an error | Align by gene name; error when panels are incompatible |
| A2 | `stability.R` skips failed refits but divides by all requested replicates | Inject a fit that fails 30% of the time; check the reported counts and frequencies | Report `n_attempted`, `n_succeeded`, failure reasons, and recovery both conditional and unconditional on success |
| A3 | `.greedy_match()` is a greedy one-to-one match | Build a 3×3 similarity matrix where the greedy and optimal assignments differ | Add `clue::solve_LSAP` with dummy rows/columns and a minimum-similarity cutoff |
| A4 | `benchmarks/coordination/methods.R` merges clusters by single linkage | Chain of three factors with A~B and B~C high but A~C low | Require every within-cluster similarity to pass the threshold (complete linkage) |
| A5 | `data-model.R` has no joint gene panel | Construct data with unequal panels; check what `$gene_universe` returns | Add a `gene_universe` and a view × gene eligibility mask (Section 4) |
| A6 | `principal-angles.R` permutes gene labels | Read the code: are labels permuted jointly or separately per view? Stratified? | Diagnostics only. Calibrate against the Section 8.3 null |
| A7 | S7 in `sim_param.R` reuses the signed loading; `same_genes_independent` reuses only the gene support | Read the code; compute the cosine between view loadings in each simulator | Label these as G1/G2-type and G4-type truths respectively in the new benchmark |

The `prereg-v1` model-selection benchmark and `benchmarks/coordination/` stay frozen under their original estimands.

---

## 4. Preprocessing contract for joint fitting (new; must be written before any fitting)

The baseline `fit_celltype_programs()` keeps its own preprocessing. The joint backend uses a separate, logged contract:

1. **Normalisation.** Library-size normalise each pseudobulk and log-transform with a fixed pseudocount (for example, log1p of CPM/100). Use the same transform in every view.
2. **Eligibility mask.** M_cj = 1 if gene j is "expressed" in view c (prespecified rule, for example ≥ 10% of pseudobulks with CPM ≥ 1 and at least the view's minimum cell count). **If M_cj = 0, the entry is missing (NA), not zero.** This keeps cell-type identity genes and genes not expressed in a view from creating false sharing or false absence.
3. **Centring within each view.** Subtract each gene's mean within its view. Without this, mean differences between cell types appear as shared factors. MOFA+ does the same for groups (MOFA FAQ).
4. **Scaling.** Default: **no per-view gene standardisation.** Scaling genes separately in each view changes the loading vector differently in each view, so an "exact core" could not exist. Residual variance is modelled instead (step 5). Sensitivity analysis: a single gene-wise scale shared by all views.
5. **Precision.**
   - Known sampling variance S²_cij ∝ (log-scale delta-method variance)/(number of cells in pseudobulk i of view c), passed to `flashier` as a fixed `S`.
   - Plus an estimated gene-wise extra variance (`var_type = 2`).
   - Limitation: `flashier` cannot estimate a separate residual variance per view and gene. If view-specific overdispersion matters, M3b adds a custom variance update.
6. **Feature selection.** Select highly variable genes on the *view-centred* data. Take the union of per-view HVG sets, capped at G_max. Selection uses only the data, never the truth.
7. **Covariates.** Regress out known technical covariates within each view before centring (MOFA FAQ recommends this). Log the design.
8. **Logging.** Write the gene universe, mask, per-view n, cell counts, transform parameters and package versions into the fit object.

---

## 5. Statistical model

### 5.1 Exact-core model, written as a stacked EBMF

Per view:

X_c = Σ_k z_ck g_kᵀ + E_c,  E_c,ij ~ N(0, S²_cij + σ²_j)

Here z_ck ∈ ℝ^{n_c} are view-specific scores and g_k ∈ ℝ^G are shared loadings. A private factor is one whose z_ck is zero in all but one view.

Stack the rows: Y = [X_1; …; X_C] ((Σ n_c) × G), l_k = [z_1k; …; z_Ck]. Then

**Y = Σ_k l_k g_kᵀ + E**, with priors

- Gene loadings: g_kj ~ G_k (point-Laplace by default, estimated by empirical Bayes; it may be point-exponential for non-negative programs).
- Stacked scores: l_k,i ~ N(0, v_k,b(i)), where b(i) is the view of row i and **v_k,b ≥ 0 is estimated by empirical Bayes**. If v_k,b = 0, factor k is exactly zero in view b.

How this covers the v1 model:
- The v1 inclusion indicator h_ck and amplitude a_ck are absorbed into v_k,c.
- Private factors are single-block factors, and shared programs are multi-block factors.
- "Zero shared programs" means every factor has one active block.
- The independent-only model is the special case where each factor lives in one block. With a single view the method reduces exactly to EBMF.
- There is **no coupling of scores across views**: rows from different views are different random variables, so score correlation cannot affect the fit. Scenario G5 therefore passes by construction (Section 7.3).

Implementation:
- Use `flashier::flash` with `ebnm_fn = c(block_normal_ebnm(blocks), ebnm_point_laplace)`.
- The custom EBNM function solves a one-dimensional maximum marginal likelihood for each block.
- The toy prototype is in `spike/block_ebnm.R`.
- `flashier` calls the function once with a length-3 probe vector, so the function must handle that input.

Identifiability:
- **Scale:** fixed by `flashier`'s LDF normalisation, with ‖g_k‖ = 1 reported.
- **Sign:** a global sign flip (g, l) → (−g, −l) is the only sign ambiguity. A per-view sign cannot be identified under a symmetric score prior, so the v1 option of "sign-constrained a_ck" is removed.
- **Rotation:** sparse non-Gaussian loading priors break rotational invariance in practice. When two programs share most of their genes, report the subspace instead (Section 9.3).
- **Duplicates:** after the backfit, merge factors whose |cos(g_k, g_k')| > 0.95 and whose active views overlap, then refit. Record every merge.

### 5.2 Membership inference: cross-fitted block test (required)

The fitted v_k,b is a point estimate on the boundary of the parameter space. Under the null it is above zero about half the time. In addition, g_k was estimated using block b's own data, so a plug-in test reuses the same data twice. The toy test confirms this: naive calls gave **65–75% false sharing on null data** (Section 6).

v2 rule (Section 6 shows its toy performance):
1. For factor k, the **anchor** is the block with the largest v_k,b.
2. For every other block b: re-estimate the loading **without block b**, f_{−b} = R_{−b}ᵀ l_{−b} / ‖l_{−b}‖², where R is the partial residual for factor k. Project block b's partial residual onto f_{−b}, using the per-gene residual precision. Compute the likelihood-ratio statistic for v_b > 0 against v_b = 0.
3. Reference distribution: a 50:50 mixture of a point mass at 0 and χ²₁ (chi-bar-square). Use a Bonferroni correction over K × (C − 1) tests, so the critical value is χ²₁ quantile at 1 − 2α/m. Default α = 0.05 (family-wise).
4. Membership A_k = {anchor} ∪ {b : test significant}. A view is **indeterminate** if a power proxy says the test could not detect the program there (for example, the expected LR at the anchor's v_k is below the critical value).

Known approximations, to be checked in M3:
- the scores l are plug-in estimates fitted with all data;
- uncertainty in F is ignored;
- Bonferroni is conservative. Later options: Holm's method, or a closed testing procedure over subsets of views.

Alternatives to benchmark in M6:
- a spike-and-slab block prior, (1 − π_k)δ_0 + π_k N(0, v_k,b), whose posterior inclusion probabilities would need their own calibration study;
- an ELBO-based nullcheck for each block (refit with v_k,b fixed at 0 and compare the ELBO).

### 5.3 Residual recruitment
In the stacked model, recruitment happens automatically: every block's data informs every factor during the greedy and backfit steps. The separate step 5 in v1 is only needed for the B1+R comparator.

### 5.4 Flexible model (M5, only if M4 passes)
w_ck = a_ck g_k + δ_ck, with δ_ck ⟂ g_k and a shrinkage prior on δ_ck. This has the same structure as LIGER's W + V_i without non-negativity. Frame it that way.

Implementation route: for each shared factor, add view-specific "deviation" factors whose scores are tied to the parent factor's scores in the same view. This requires a custom update, because `flashier` cannot tie scores across factors. Start with fixed shrinkage, and add an empirical Bayes prior on δ only if fixed shrinkage helps.

### 5.5 Hierarchical extensions (M6)
- Share the prior on v_k,· across factors (a hierarchical prior on block variances).
- Spike-and-slab membership.
- View-specific gene-wise residual variance.

Each one goes in only if it beats M5 on the preregistered endpoints.

---

## 6. Feasibility evidence (toy test, run 2026-10-10)

**Setup:**
- Tier-1 Gaussian data: C = 4 views, n = 60 per view, G = 1000 genes, noise σ = 1.
- One private sparse factor per view: 60 active genes, r = 3.
- One shared sparse signed core in views 1–3, with independent scores in each view.
- **r** = spike variance divided by the single-view BBP threshold √(G/n).
- Methods:
  - B0: independent `flashier` per view (normal scores, point-Laplace loadings, `var_type = 2`, backfit).
  - B1+R: each B0 loading tested in the other views with the same Bonferroni-corrected LR test. The candidate closest to the truth was chosen by an oracle, which favours B1+R.
  - Stacked: block-normal score prior, point-Laplace loadings.
- Code and per-replicate results: `spike/` (`block_ebnm.R`, `fit_compare.R`, `spike_results*.csv`).

**Run 1: naive membership calls (20 replicates per row)**

| Scenario | Stacked, v_b > 0 in ≥ 2 views | Stacked, plug-in LR | B1+R (uncorrected) |
|---|---|---|---|
| G0 null | 0.75 | 0.45 | 0.15 |

**Run 2: cross-fitted Bonferroni membership (40 null replicates, 20 for the others)**

| Scenario | r | Stacked datasets with a sharing call (95% CI) | B1+R with a sharing call (95% CI) | Stacked cosine to g | B0 cosine to g | Stacked exact membership | B1+R exact membership | Private cosine, stacked / B0 |
|---|---|---|---|---|---|---|---|---|
| G0 null | – | **0.00 [0.00, 0.09]** | 0.03 [0.00, 0.13] | – | – | – | – | 0.79 / 0.85 |
| G2 partial | 0.5 | 0.15 [0.03, 0.38] | 0.10 [0.01, 0.32] | 0.09 | 0.03 | 0.00 | 0.00 | 0.82 / 0.87 |
| G2 partial | 1.0 | 0.90 [0.68, 0.99] | 0.10 [0.01, 0.32] | **0.73** | 0.05 | **0.80** | 0.10 | 0.83 / 0.86 |
| G2 partial | 2.0 | 1.00 [0.83, 1.00] | 0.95 [0.75, 1.00] | 0.94 | 0.47 | 1.00 | 0.95 | 0.86 / 0.87 |
| G4 same genes, different weights | 2.0 | 0.15 [0.03, 0.38] | 0.10 [0.01, 0.32] | 0.10 | 0.07 | – | – | 0.80 / 0.86 |

**What this shows (toy only):**
1. **The stacked formulation runs in `flashier` 1.0.7 without modification.** One fit takes about 1–2 s at this size.
2. **The borrowing window is real and where theory predicts.**
   - At r = 1, a single view is at its detection threshold. B0 finds nothing, and B1+R has no candidate to recruit.
   - The stacked fit pools 3 views (threshold factor 1/√3 ≈ 0.58), recovers the core, and gets exact membership in 80% of datasets.
   - At r = 2, both approaches work.
   - At r = 0.5, which is below even the pooled threshold, neither works.
3. **Naive membership is badly miscalibrated.** The cross-fitted test is required (Section 5.2).
4. **Two warning signs need checking at scale:**
   - (a) G4 (same genes, unrelated weights) gave 15% sharing calls against 0% under G0. The CIs overlap, but this may be leakage driven by shared gene support.
   - (b) Private-factor recovery is consistently a little lower for the stacked fit (by about 0.01–0.06 in cosine), so H2 is not guaranteed.

**Not shown:** realistic noise, heteroscedastic precision, unequal n, more views, technical signatures, count data, or real data. All of these are left to M1–M7.

---

## 7. Simulation design

### 7.1 Tiers
- **Tier 1, Gaussian parametric** (Section 7.3). Used to debug, calibrate and test H1–H3.
- **Tier 2, count-based.** scDesign3 (Song et al. 2023) or Splatter fitted to real PBMC data. Vary cells per donor, gene mean and dispersion, and cell-type prevalence. Pseudobulk and preprocess through the Section 4 pipeline.
- **Tier 3, semi-synthetic.** Inject programs into COMBAT/Stephenson pseudobulk. This is realistic, but the background is not a null.
- **Tier 3N, permutation null on real data** (Section 8.3).

### 7.2 Signal axis
For a shared program with membership m = |A| and per-view n:
- single-view threshold τ₁ = √(G/n);
- pooled threshold τ_m ≈ √(G/(m·n)).

Set r = θ/τ₁ ∈ {0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0}. Expected regimes:
- **r < 1/√m:** below both thresholds (sparsity priors may lower the effective threshold);
- **1/√m ≤ r ≤ ~1.5:** the **borrowing window**, where the primary endpoint is evaluated;
- **r > ~1.5:** both approaches should work.

The window boundaries are computed from the design in each cell, not tuned afterwards.

### 7.3 Scenario families (Tier 1)

| ID | Gene-loading truth | Scores across views | Expected outcome | Notes |
|---|---|---|---|---|
| G0 | Private factors only | Independent | No sharing | Overlapping sparse gene supports and correlated gene background (hard null) |
| G1 | One exact core in all views | Independent | Fully shared | |
| G2 | Exact core in a subset | Independent | Correct partial membership | Main H1 scenario |
| G3 | Core plus controlled deviation (cos ∈ {0.95, 0.8, 0.5}) | Independent | Shared down to some cosine; this sets the limits of the exact-core model | Exact-core performance should fall gracefully |
| G4 | Same gene support, independent signed weights | Independent | No sharing call | Toy warning sign (a) |
| G5 | Different programs | Coordinated (matched donors, ρ ∈ {0.5, 0.9}) | No sharing call | Passes by construction for stacked/B1. Test it to confirm, and so that MOFAcellulaR contrasts are fair |
| G6 | Correlated but distinct loadings (cos 0.3–0.7) | Independent | No merging | |
| G7 | Shared technical signature (library size, stress) on top of private biology | Either | Recurrence detected, then flagged as nuisance | Counted separately from false sharing |
| G8 | One program split across two factors, plus overlapping programs | Independent | Recovery that accounts for splits and merges | |
| G9 | Shared core, weak in one member view (r_view ≈ 0.5) | Independent | Weak view recruited without forcing others in | |
| **G10** (new) | Mean differences between cell types left uncentred, plus identity genes | Independent | Tests the preprocessing contract | Run once with and once without step 4.3 |
| **G11** (new) | Heteroscedastic precision (cells per donor vary 10–500) | Independent | Calibration holds | Uses `S` |
| **G12** (new) | Unequal n_c (20–300) and missing views | Independent | Calibration holds; indeterminate calls reported | |

### 7.4 Factor grid
- **Pilot (P0):**
  - C ∈ {4, 8}; n = 60; G = 2000.
  - 1 shared program; 2 private programs per view; 60 active genes.
  - r over the Section 7.2 grid; scenarios G0, G2, G4, G9.
  - 50 replicates per cell.
- **Confirmatory:**
  - C ∈ {4, 8, 16}; n ∈ {30, 100}; G ∈ {2000, 5000}; 0–3 shared programs; 1–3 private programs per view; 30–150 active genes.
  - All scenarios G0–G12.
  - 200 replicates for positive cells and 300 for null cells (Section 9.4).
  - **New seeds, never used in the pilot.**
- **Stress:** C = 20; n = 20–300; G = 10,000; 70% missing views.

### 7.5 Self-checks before any method runs
Automated assertions on every simulated dataset:
- membership sets match the specification;
- realised loading cosines and support overlap match the targets;
- realised score correlation between views |ρ| < 0.1 in independent designs;
- realised r is within 10% of the target;
- masks and missing views behave correctly;
- no unplanned shared structure in nulls (top singular vectors of the stacked null match the per-view spectra);
- results reproduce under a fixed seed.

---

## 8. Methods and comparators

### 8.1 Methods

| Label | Method | Implementation notes |
|---|---|---|
| B0 | Independent `flashier` per view | `fit_celltype_programs()` defaults, unchanged |
| B1 | B0 + signed consensus clustering of loadings | Optimal matching, complete linkage, coherence screen |
| **B1+R** | B1 + residual recruitment of each consensus candidate into every view, with the same Bonferroni-corrected test | **Main ablation:** separates the gain from joint fitting from the gain from recruitment |
| B2 | LIGER iNMF / UINMF | Non-negative input (shifted log or counts). Shared W versus specific V_i compared by cosine |
| **B3** | **MOFA+ multi-group**: one view, cell types as groups | Group-wise ARD gives per-group factor activity. Membership = groups with explained variance above a threshold calibrated on G0 |
| B3b | MSFA / Tetris | Small grid only (C ≤ 4, G ≤ 500 after HVG selection), because MCMC is slow |
| **B4** | Stacked EBMF with block-normal score prior and cross-fitted membership | Section 5.1–5.2 |
| B5 | B4 + view-specific deviations, fixed shrinkage | M5 |
| B6 | B5 + hierarchical empirical Bayes | M6 |
| Oracle | True memberships given; loadings estimated | Upper reference |

MOFAcellulaR is **not a gene-sharing comparator**. It is used only in the post-hoc coordination analysis on matched donors.

### 8.2 Fairness
- Same input matrices after the Section 4 pipeline. Method-specific transforms are logged (for example, non-negative input for LIGER).
- Same rank budget: K_max = 2 × true total rank, or ≥ 20.
- Same tuning budget, with tuning done on pilot seeds only.
- Matched-donor and unpaired tracks are run separately.
- Each method's membership threshold is **calibrated on G0 to the same null FWER** before power is compared (Section 9.2).

### 8.3 Real-data null and controls
- **Stratified within-view gene permutation:**
  - In each view separately, permute gene labels within bins of mean × variance (20 × 5 bins).
  - Within each view this only relabels genes, so the spectrum, sparsity and noise are unchanged.
  - Correspondence between views is destroyed, so any sharing call is false, apart from artefacts linked to the mean–variance bins.
  - Use it to calibrate the real-data FWER.
- **Positive controls (expected shared):**
  - interferon response (ISG15, IFI6, MX1, …);
  - cell cycle (G2M/S genes in proliferating compartments);
  - sex genes (XIST, RPS4Y1, DDX3Y). In mixed-sex cohorts these act as a shared program whose scores are also coordinated;
  - dissociation/stress (FOS, JUN, HSPA1A). This is shared but technical, so it is flagged.
- **Negative controls (expected private):** cell-type identity programs restricted to one view (for example, a B-cell immunoglobulin program in non-B views should be masked or private).

---

## 9. Endpoints, metrics and statistics

### 9.1 Matching
- Hungarian assignment (`clue::solve_LSAP`) on |cos| between true and estimated programs, with dummy rows/columns and a minimum |cos| of 0.5 to count as a match.
- **Split-aware recovery:** ‖P_Ŝ g_k‖, where P_Ŝ is the projection onto the span of estimated loadings with |cos| > 0.3 to g_k. Report this next to the 1:1 cosine.

### 9.2 Endpoints

| Role | Endpoint | Definition |
|---|---|---|
| **Primary (H1)** | Exact-membership recall in the borrowing window | Fraction of true shared programs that are matched (cos ≥ 0.8) **and** have Â_k = A_k. Measured in G2/G9 cells with 1/√m ≤ r ≤ 1.5, at each method's G0-calibrated 5% FWER threshold. Hypothesis: B4 − (B1+R) ≥ **0.10**, with the lower 95% CI bound > 0 (paired by replicate) |
| **Co-primary (H2)** | Private-program preservation | Mean matched cosine for private programs. Non-inferiority: B4 − B0 ≥ **−0.03** (lower 95% CI) |
| **Co-primary (H3)** | Null FWER | Fraction of G0 (and Tier-3N) datasets with any sharing call. Upper 95% CI ≤ 0.08 at the default threshold |
| Secondary | Membership F1; matched loading cosine; shared-program recall versus r curves; gene-support AUPRC; split/merge counts; G4/G6 false-merge rate; G7 nuisance flag accuracy | |
| Secondary | Power at matched FWER (ROC-style curve of membership recall against G0 FWER as α varies) | Separates power from calibration |
| Secondary | Held-out prediction: masked entries reconstructed for held-out donors; projection from source to target views | Supporting only. Better reconstruction does not mean sharing was called correctly |
| Engineering | Runtime, peak RAM, convergence failures, seed sensitivity (adjusted Rand index of memberships across 5 initialisations) | |

The ±0.10 and −0.03 margins are **proposed values, to be confirmed after the pilot** and before the preregistration is frozen.

### 9.3 Ambiguity reporting
When two estimated programs have |cos| > 0.7 and their memberships overlap, report them as a shared subspace with its principal angles, not as two separate programs.

### 9.4 Monte Carlo sizes
- **Null FWER:** at a true rate of 0.05, 300 replicates give a 95% CI half-width of about ±0.025. If 0/300 are observed, the upper bound is about 0.012.
- **Paired recall difference:** a difference of 0.10 with about 25% discordant pairs needs about 200 replicates for 80% power (McNemar test, two-sided α = 0.05).
- **Failures:** report attempted, succeeded and failed fits. A failed fit counts as a **miss** for recall and is excluded from FWER, with the failure rate reported. It never counts as a true negative.

### 9.5 Leakage control
- Rank, priors, thresholds and matching are chosen on training donors or pilot seeds only.
- Held-out evaluation uses donor-level splits.
- Cross-view transfer: learn g_k in the source views, fit target scores on target training genes, and evaluate on masked target genes.

---

## 10. Real-data plan (M7)

1. **COMBAT** (discovery) and **Stephenson et al.** (replication), PBMC from COVID-19 cohorts. 8–12 major cell types, donor-level pseudobulks.
2. Run B0, B1+R, B3 and B4 with the same preprocessing.
3. **Replication:** for each shared program in COMBAT, find the best |cos| in Stephenson and compare the memberships (Jaccard index), against the stratified-permutation null distribution.
4. **Controls** as in Section 8.3. Expected: interferon shared broadly across myeloid and lymphoid views; identity programs private; stress programs flagged.
5. **Post-hoc analyses:**
   - donor coordination of the scores of each shared program across views (where donors are matched);
   - association with disease severity using the per-view scores;
   - comparison with MOFAcellulaR multicellular factors.
6. **Scale-up:** HLCA (tissue, many cell types, unequal n) and OneK1K (about 1k donors, test of n).
7. **Small cell types:** check that indeterminate calls appear where expected, rather than false absences.

Gate: positive controls recovered, the stratified null is calibrated, programs replicate across cohorts above the null, and no unexplained technical programs remain.

---

## 11. Milestones and gates (revised order)

| M | Content | Gate | Rough effort |
|---|---|---|---|
| **M0** | Audit tests A1–A7 (Section 3); fixes; `gene_universe` and mask; Section 4 preprocessing contract with logging; provenance (versions, seeds) | Existing tests pass; baseline numbers unchanged; A1–A7 tests pass | 1–2 wk |
| **M1** | Tier-1 simulator for G0–G12 with Section 7.5 self-checks; signal parameterised by r | All invariants pass before any method is run | 1–2 wk |
| **M2** | Package the block-prior EBNM function (`R/block-ebnm.R`, with unit tests: probe-call handling, v = 0 boundary, recovery of known v); stacked fit wrapper; cross-fitted membership; B1 and B1+R | Unit tests pass; toy results reproduced | 1–2 wk |
| **M3** | Calibration study on G0/G4/G11/G12; refine the membership test (Holm, closed testing, power proxy for indeterminate calls); duplicate merging; multiple initialisations | FWER ≤ 5% on G0 **and** G4 (upper CI ≤ 8%); seed ARI ≥ 0.9 | 2–3 wk |
| **P0** | **Pilot** (Section 7.4) with B0, B1, B1+R, B3 (MOFA+), B4; confirm the margins; write `benchmarks/gene_sharing/PREREG.md`; tag `prereg-gs-v1` | Preregistration frozen | 2 wk |
| **M4** | Confirmatory Tier 1 (all methods except B5/B6) plus Tier 2 subset; primary analysis | **H1, H2 and H3 all pass, otherwise STOP** (write up the negative or partial result) | 3–4 wk (HPC) |
| M5 | Flexible deviations (Section 5.4), fixed shrinkage | Beats B4 on G3/G6/G9 with no loss on G0 | 3–4 wk |
| M6 | Hierarchical empirical Bayes options (Section 5.5) | Beats M5 on preregistered secondary endpoints | 4+ wk |
| M7 | Real data (Section 10) | Section 10 gate | 3–4 wk |
| M8 | API, documentation, vignettes, benchmark report | `R CMD check` clean; vignette reproduces the main figures | 2 wk |

---

## 12. Software design

### 12.1 Files
```text
R/
  block-ebnm.R              # block_normal_ebnm(), print/summary for fitted g
  gene-universe.R           # build_gene_universe(), eligibility mask, NA semantics
  shared-preprocess.R       # Section 4 contract, logged
  shared-fit.R              # fit_shared_programs(): stacking, flash call, merging
  shared-membership.R       # crossfit_membership(), power proxy, multiplicity
  shared-extract.R          # shared_programs(), program_membership(), private_programs()
  consensus-programs.R      # B1 / B1+R
  matching.R                # optimal assignment, split-aware recovery
benchmarks/gene_sharing/
  PREREG.md
  sim/ parametric.R  count_based.R  semisynthetic.R  permutation_null.R
  methods/ b0.R b1.R b1r.R liger.R mofa_multigroup.R msfa.R b4.R
  evaluate.R  summarize.R  report.Rmd
tests/testthat/
  test-block-ebnm.R  test-gene-universe.R  test-shared-fit.R
  test-membership-null.R  test-membership-recovery.R  test-matching.R
  test-sim-gene-sharing.R  test-adversarial.R
```
"Hierarchical" is not used in any user-facing name until M6 delivers a hierarchical model.

### 12.2 API
```r
pp  <- prepare_shared_input(data, gene_rule = "pct10_cpm1", center = "view",
                            scale = "none", precision = "cells")
fit <- fit_shared_programs(pp, max_factors = 30, loading_prior = "point_laplace",
                           membership = list(method = "crossfit", alpha = 0.05,
                                             correction = "bonferroni"),
                           n_init = 5, seed = 1)
shared_programs(fit)      # g_k, membership, class, ambiguity flags
program_membership(fit)   # program x view: supported / unsupported / indeterminate / unobserved
private_programs(fit)
validate_sharing(fit, null = "stratified_permutation", n_perm = 100)
```
The fit object stores: g_k (posterior mean, SD, lfsr), per-view scores, v_k,b, membership statistics and decisions, merges, ELBO trace, preprocessing log, and versions/seeds.

### 12.3 Adversarial tests (unit-level)
- Null data initialised with shared candidates → no sharing.
- A view with only noise → no membership in that view.
- Duplicated cores → merged.
- Global sign flip → identical programs.
- Permuted view order → identical up to labels.
- Reordered genes and unequal panels → aligned or a clear error.
- A gene missing in one view → NA, never treated as zero.
- A weak program in one view → recruited.
- Coordinated scores with different genes → no sharing.
- Technical signature → shared and flagged.
- `flash` failure → reported as failure, never as "no sharing".

---

## 13. Preregistration contents (`benchmarks/gene_sharing/PREREG.md`)
1. Estimand and membership definition (Section 1.3).
2. Preprocessing contract and versions (Section 4).
3. Methods, settings and tuning budgets (Section 8).
4. Scenario grid, seeds and replicate counts (Section 7.4, Section 9.4).
5. Primary and co-primary endpoints, margins and the decision rule (Section 9.2).
6. Matching and failure handling (Section 9.1, Section 9.4).
7. Secondary and exploratory analyses, labelled as such.
8. Stop/go rule at M4.
9. Deviations log.

---

## 14. Risks and mitigations

| Risk | Signal | Mitigation |
|---|---|---|
| Membership test anti-conservative at realistic scale (plug-in scores, gene correlation) | G0/G11/Tier-3N FWER > 5% | Cross-fit scores as well; closed testing; calibrate on stratified permutations; report empirical rather than nominal α |
| Shared gene support leaks into sharing calls (toy G4 = 0.15) | G4 FWER above G0 | Higher cosine cutoff for matching; test G4 explicitly in M3; consider a test on signed loading agreement |
| Private-factor recovery degrades (toy −0.01 to −0.06) | H2 fails | Higher K_max; backfit from B0 initialisation; refit private factors per view after joint fit |
| Per-gene residual variance shared across views is misspecified | Poor calibration in G11 | Fixed `S` from cell counts; custom view × gene variance update (M3b) |
| MOFA+ multi-group already achieves H1 | B3 ≈ B4 | The contribution becomes calibrated membership plus sparse empirical Bayes loadings; report honestly |
| Tetris outperforms B4 on small problems | B3b > B4 at C ≤ 4 | Position B4 as the scalable alternative; quantify the runtime trade-off |
| Real-data programs are mostly technical | Stress or ribosomal programs dominate | Covariate regression, nuisance flags, sensitivity analysis excluding technical gene sets |

---

## 15. Decisions needed from you

1. **Primary-endpoint margins:** is +0.10 exact-membership recall over B1+R and −0.03 private non-inferiority acceptable as the starting values for the pilot?
2. **Scaling default:** no per-view gene scaling (recommended), or a single gene-wise scale shared by all views?
3. **Gene-eligibility rule:** is ≥ 10% of pseudobulks with CPM ≥ 1 right for your cohorts?
4. **Competitor scope:** run Tetris/MSFA on the small grid, or cite them only?
5. **Matched-donor track:** do you want the G5 coordination contrast against MOFAcellulaR in this workstream, or keep it in the existing coordination benchmark?

---

### Key references
Argelaguet et al. 2020 *Genome Biol* (MOFA+). De Vito et al. 2021 *Ann Appl Stat* (BMSFA). Grabski et al. 2023 *Ann Appl Stat* (Tetris). Gaynanova & Li 2019 *Biometrics* (SLIDE). Lock et al. 2013 *Ann Appl Stat* (JIVE). Kriebel & Welch 2022 *Nat Commun* (UINMF). Kunes et al. 2023 *Nat Biotechnol* (Spectra). Kotliar et al. 2019 *eLife* (cNMF). Andreatta & Carmona 2024 (GeneNMF). Ramirez Flores et al. 2023 *eLife* (MOFAcellulaR). Jerby-Arnon & Regev 2022 *Nat Biotechnol* (DIALOGUE). Song et al. 2023 *Nat Biotechnol* (scDesign3). Baik, Ben Arous & Péché 2005 *Ann Probab* (BBP phase transition). Wang & Stephens 2021 *JMLR* (EBMF); `flashier` 1.0.7 documentation.
