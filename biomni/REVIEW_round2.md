# Review round 2 (revised after code review): replication, corrected metrics, principal angles vs PMD

Supersedes the first version of this file. All code and outputs: `benchmarks/biomni_replication/`
(repo package untouched). Pseudobulk from CELLxGENE Census 2025-11-08. Fork = Biomni's `cellprograms/`
plus `integrate-pmd.R` (patched only to route `method = "pmd"`). Nothing below is committed.

## What changed since the first version
A code review of my own work found three defects that invalidated earlier conclusions. All are fixed and everything was re-run.
1. **Wrong input matrix.** I had pseudobulked `raw.X` (raw counts). Biomni's trace is only reproducible from the
   log-normalized `X`: sums of `X`, truncated to integers (`floor`), then Biomni's CPM/log2/top-3000 recipe. The `floor` is load-bearing for replication: without it
   marker values differ by up to 0.84 and only 49/120 marker cells match within 0.05. Keep it when replicating; it is not a recommended default for new analyses.
   With that input the COMBAT pipeline reproduces Biomni almost exactly (section 1). The "irreducible ±0.05 baseline gap"
   I reported earlier was caused by this and by z-scoring, not by missing information.
2. **`cv_ridge` bug** (leverage omitted the intercept term, so high-dimensional representations got a near-zero penalty
   and R2 of -22 to -85). Fixed; test added. It had hurt cp specifically.
3. **Hand-counted verdict table.** Replaced by `summarize_eval.R` -> `results/eval/verdicts.csv` (reproducible).
Also corrected: PA-vs-PMD wording (section 3), a mislabelled scenario row, and hygiene (.gitignore, paths, env note).
Old raw-count outputs are kept under `data/*_raw`, `results/fits_rawcounts`, `results/eval_rawcounts`, `results/real_rawcounts`.

## 1. Replication of Biomni on COMBAT
| check | result |
|---|---|
| marker table (trace cell 4, 12 markers x 10 cell types, mean log2CPM) | all 120 cells within 0.05 (max abs diff 0.00) |
| top-3000 gene panels vs Biomni's program top genes | 100% overlap in all 10 cell types (raw-count arm: 59-95%) |
| donor counts per cell type | exact for COMBAT, Stephenson, HLCA; OneK1K 3/9 exact (cap draw not recoverable; 1-9 donors off) |
| baselines, continuous (max abs Spearman) | composition, global/perct/grouped PCA exact (diff 0.000); clr off 0.033 (pseudocount unknown); mean 0.007 |
| baselines, categorical 5-NN macro-F1 (no z-scoring) | mean abs diff 0.022 (PCA baselines match to 3 decimals on Death28, sex, Institute; composition 0.014; clr 0.043) |
| cp_L0, Biomni config (K=30, fork defaults) vs Biomni's reported cp_L0 | mean abs diff 0.004 over 9 covariates (7 exact; worst Source -0.025) |
Biomni did not z-score the baseline representations; my first implementation did (0.064 error). Default is now unscaled.
Still unreplicated: Stephenson/HLCA/OneK1K cp numbers (Biomni used K=25 and `prior_max_iter=3` there; I ran K=30), `random` (unknown seed), clr pseudocount.
Stephenson Age is not in the source h5ad; I derived it from `development_stage` (NaN for 62/120 donors), so it is not comparable.

## 2. Corrected benchmark (`verdicts.csv`)
Method: categorical = available-case kNN macro-F1 minus the masked dimension-matched random floor; continuous =
max |Spearman| minus a floor at the representation's own dimension; best baseline chosen by excess (not raw value); tie band ±0.05;
a second verdict against the median baseline removes the selection bias of taking the max over 4-5 baselines.
Excluded: HLCA `age_or_mean_of_age_range` (duplicate), Stephenson `Collection_Day`/`Resample` (degenerate).

Tally over biology covariates, summed over four cohorts, cp / tie / baseline (vs best baseline | vs median baseline):
| config | covariate class | n | vs best | vs median |
|---|---|---|---|---|
| K=30 (Biomni) | expression-driven (sex, Age, BMI, TimeSinceOnset) | 10 | 2 / 4 / 4 | 4 / 4 / 2 |
| K=10 (repo default) | expression-driven | 10 | 3 / 3 / 4 | 7 / 2 / 1 |
| K=30 | compositional (Source, Outcome, Death28, Status, Smoker) | 10 | 0 / 1 / 9 | 0 / 1 / 9 |
| K=10 | compositional | 10 | 0 / 3 / 7 | 1 / 7 / 2 |

Reading:
- **K=10 beats Biomni's K=30** on COMBAT sex (excess 0.51 vs 0.33), Source (0.31 vs 0.18), Age (0.23 vs 0.16). The repo default is the right config.
- **Sex is the robust cp win** (COMBAT K=10 0.51 vs best baseline 0.32; Stephenson 0.49 vs 0.29; HLCA 0.39 vs 0.11; OneK1K ties).
- **COMBAT Age:** cp beats the dimension-matched floor by 0.23 (K=10) vs best baseline perct_pca 0.21: a tie vs the best baseline and a win vs
  the median. Cross-validated ridge R2 above its floor (fixed): cp 0.39 (K=10) / 0.28 (K=30) vs best baseline 0.16. So the first-version doubt about COMBAT Age only holds for the
  K=30 max-Spearman comparison, and mostly dissolves under ridge.
- **Compositional covariates go to composition / clr / PCA baselines**, as Biomni also found. Stephenson Age and OneK1K age go to composition/clr.
- **Zero-fill confound is real but smaller than first reported.** COMBAT Institute kNN, K=10: 0.677 zero-filled vs 0.532 available-case
  (floor 0.473); perct_pca shows the same drop (0.673 -> 0.532). At K=30 cp sits at the floor (0.470). So some Institute information remains in
  K=10 cp after removing missingness (+0.06 over floor).
- **Panel B** (excess over floor on technical covariates, lower is better): cp is cleaner than the most confounded baseline on 8 of 10 rows, and level
  on 2; against the median baseline it is cleaner on 3, level on 6, worse on 1 (Stephenson Site, K=10).
- Not done: block integration at K=10, GloScope/PILOT, stability refits, cross-cohort matching.

## 3. Principal angles vs PMD
Simulation: 20 seeds, S1-S9, missing 0/0.3 (`results/pa_vs_pmd_sim_summary_v2.csv`). S9 = S3 plus a donor-level nuisance acting through genes in every cell type.
Caveats: the PA and PMD "shared" rules are not matched in error rate (PA: any of up to min(k_a,k_b) directions, no multiplicity correction; PMD: one test per MCP then r > 0.5),
so TPR/FPR comparisons mostly reflect the rules; AUROC is the fairer comparison. S2/S4/S6 truth shares every pair, so their TPR = 1 says little. S5 (r = 0.7) and S7 (shared genes)
"false positives" depend on the truth definition. Nuisance fixes below use the true nuisance (oracle) unless stated.

| | PA-scores | PA-loadings | PMD |
|---|---|---|---|
| AUROC S3 (one shared pair), missing 0 / 0.3 | 1.00 / 1.00 | 0.68 / 0.68 | 1.00 / 1.00 (continuous score; 0.85 at 0.3 when non-significant pairs are zeroed) |
| AUROC S9 raw, 0 / 0.3 | 0.86 / 0.82 | 0.66 / 0.67 | 0.70 / 0.64 |
| AUROC S9 after oracle residualization | 1.00 / 1.00 | n/a | 1.00 / 1.00 (0.95 / 0.80 with significance gating) |
| AUROC S9 with ESTIMATED nuisance | 1.00 / 0.93 | n/a | 0.81 / 0.77 |
| FPR all-private S1/S7 | 0.04-0.09 | 0.05-0.12 | 0.00 |
| FPR S9 raw | 1.00 | 0.07 | 1.00 |
| FPR S9 oracle-residualized | 0.02-0.05 | n/a | 0 (TPR 0.9 / 0.6) |
| FPR S9 estimated nuisance | 1.00 | n/a | 1.00 |
Estimated nuisance = first axis shared across cell types; correlates 0.88 / 0.82 with the true nuisance on S9 and only 0.68 on S8.
Reading: PA ranks pairs equal or better than PMD in every scenario. PMD is more conservative at its default call rule. Both call every pair shared when a
nuisance acts in every cell type, and residualization with the exact nuisance fixes both; with an estimated nuisance the ranking mostly survives (PA) but the calls do not (FPR 1.0).
Residualizing fitted scores, not data, leaves stage-1 capacity spent on the nuisance.

Real COMBAT, K=10, Biomni-faithful input (`results/real/`): PA calls 40 of the 40 callable pairs shared (40/40 after residualizing on sex + Institute); 5 of the 45 pairs cannot be called
(DC-GDT, DC-PB, DP-GDT, DP-PB, GDT-PB: 18-38 shared donors, reference cosine near or above 1). PMD calls 45/45 both ways. All 5 MCPs sit at the permutation
floor (1/101); with a value-shuffle null that is the resolution limit, and it means the test rejects everything (the null does not
separate coordination from any shared donor-level variable), not that there is no signal. Neither tool discriminates on real data.
With this input PMD reproduces Biomni's five MCPs: Ig V-genes, sex (Y genes), interferon, immediate-early/JUN (Institute p 6e-4, Source p 7e-17), GSTM1/IGHV.
After residualizing on sex + Institute:
- sex association disappears (p >= 0.69), though Y genes still lead some traced loadings;
- the immediate-early MCP is gone, so with this input it is consistent with being Institute-driven (its Source p 7e-17 is confounded with Institute; not tested further);
- interferon survives; the Ig/GSTM1 MCPs survive (no sex, Institute or Source association; ambient RNA or germline variation is a guess, untested);
- a new MCP appears with strong Source association (Source p 2e-11, Outcome p 2e-10; CLDN5/SOCS/IFNGR2 in CD4).
Note: on the earlier raw-count input the immediate-early axis survived residualization. PMD conclusions depend on the input normalization.

## 4. Port decision (unchanged in direction, sharpened)
- Port now: `principal-angles.R`, `stability.R`, `extract.R`, with (a) the "shared" rule `sqrt(k_a)+sqrt(k_b))/sqrt(n)` replaced by an explicit "uncallable" flag
  (blocks 5 of 45 pairs on COMBAT), (b) covariate residualization as a documented step, (c) multiplicity control across pairs.
- PMD: do not port yet. PA ranks as well or better, PMD's permutation test rejects everything (uncalibrated null), and its axes need residualization to be interpretable.
- Block integration, evbidifac, L1/L2 borrowing: unchanged; block needs a K=10 re-run first.
- Use K=10 `point_laplace`, `var_type = 1` (repo defaults); Biomni's K=30 understates cp.

## 5. Limits
- OneK1K donor cap and HLCA/OneK1K label choices are inferred (`README.md` in the work folder); Stephenson Age is reconstructed.
- cp numbers for Stephenson/HLCA/OneK1K were not matched to Biomni's (different K and prior iterations).
- clr pseudocount and `random` seed unidentified; ridge verdicts are secondary (`cv_*` columns only).
- Simulation: one simulator, 20 seeds, AUROC over 10 pairs; PA/PMD call rules not error-rate matched.
- Genotype interpretation of the Ig/GSTM1 MCPs is untested.
- Queue note: the `medium` partition QOS caps a user at 4 CPU / 8G, so jobs were moved to `all`.

## 6. External consult (Biomni, task tsk_0116fTRqFPmWW8S79CdeiPXE) and how it was used
Biomni's ranked suggestions: covariate residualization before EBMF; a calibrated sharing test (per-cell-type donor permutation, BH over pairs, explicit "underpowered" output);
missingness-aware representations; K chosen by permutation-calibrated variance and stability pruning; CV ridge-logistic and per-program association tests as default metrics;
reference mapping, case-control mixed models and cross-cohort matching with a null as extensions; MOFAcell, scITD, DIALOGUE and mc-ASTRA as comparators.
Accepted corrections (folded in above): 40/40 callable pairs, floor-p means rejects-everything, Institute wording.
Not accepted: "sex rests on one cohort" (sex won in COMBAT, Stephenson and HLCA; OneK1K tied) and "drop the integer floor" (it is needed to replicate).
Note on missingness: our per-cell-type fits already use only observed donors; zero-fill occurs in the concatenated representation and in downstream tests, so the fix is an available-case representation, not a change to the factorization.
