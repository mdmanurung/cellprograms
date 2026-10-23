# DIVAS COVID-19 Case Study: cellprograms Reanalysis

**Date**: 2026-10-13 · **Package**: cellprograms v0.2.0 · **PR**: #5

## What was done

Full reanalysis of the DIVAS COVID-19 multi-omics case study
(github.com/ByronSyun/DIVAS_COVID19_CaseStudy) with `cellprograms`, and a
systematic program-level comparison against the precomputed DIVAS results
(Zenodo record 17430294).

**Cohort**: 114 samples (85 patients: 29 paired T1/T2, 56 single-timepoint);
40 mild / 40 moderate / 34 severe; outcome = WHO ordinal score (wos).
**Blocks**: 4 scRNA-seq cell types (CD4 T, CD8 T, CD14 monocytes, NK;
8,634 shared genes) + proteomics (481) + metabolomics (763).

**Fitting**: cellprograms at DIVAS parity (center=TRUE, scale=FALSE, no
covariates, max_factors=30, point-Laplace prior, backfit+nullcheck, seed=1)
→ K=30 programs per block (180 scores). DIVAS reference: 91 components
(44 individual / 46 partial / 1 joint).

## Headline results

| Aspect | Result |
|---|---|
| Concordance | 23/91 DIVAS components matched at \|ρ\|>0.5 (51/180 reverse); 42 mutual-best pairs (19 >0.5); all 120 gene-block programs shared-classified, best-matching DIVAS components span 96 partial / 18 individual / 6 joint |
| Breadth trend | Median best-match correlation rises 0.38 (1 block) → 0.80 (5 blocks); joint component 0.78 |
| Clinical (wos) | CV ridge R²: cellprograms **0.64** vs DIVAS 0.42 (severity-only baseline 0.90); BH<0.05 associations 44/180 vs 8/91 |
| Effect-size agreement | ρ=0.17 over 91 pairs overall; significant pairs sit mostly above the diagonal (cellprograms counterpart carries ≥ as much wos signal). Top-100 loading Jaccard for 40 mutual-best gene-block pairs: median 0.16 (IQR 0.14–0.22) — scores agree, gene-level loadings differ |
| Coordination | 148 edges → 6 CMPs (marginal); 123 edges → 6 CMPs (within-severity) |
| Sharing | 120 gene-block programs in 31 groups; all 120 classified shared (mashr, frac_top_shared ≥ 0.98) |
| Pairing robustness | T1-only: 6/6 CMPs retained (50/148 edges); patient-blocked permutation: 121/148 edges (82%) survive at p<0.05 |

## Interpretation

1. **Cross-block structure agrees.** DIVAS's partial/joint components — its
   headline multi-omics finding — map onto cellprograms' shared programs
   with increasing fidelity as breadth grows. The two pipelines recover the
   same dominant biology (a severity-driven axis present across cell types
   and omics layers).
2. **cellprograms adds per-cell-type resolution.** 180 block-specific scores
   vs 91 pooled components; cellprograms scores predict wos ~50% better and
   detect 5× more significant associations, while 13 programs have no DIVAS
   counterpart (block-specific structure DIVAS pools away).
3. **Sharing is more permissive than DIVAS taxonomy.** All 120 gene-block
   programs classify shared vs DIVAS's 44/46/1 individual/partial/joint
   split — expected here because the 4 immune cell types are strongly
   correlated in this cohort and mash's lfsr is liberal under correlated
   conditions. Documented as a caveat, not hidden.
4. **Coordination is robust to the pairing structure.** 6/6 CMPs persist at
   T1-only; 82% of edges survive a patient-blocked permutation null.

## Files

- `divas_case_study/` — 7 figures (PNG) + 12 result tables (CSV) + REPORT.md
- Repo: `benchmarks/divas_comparison/` (scripts 01–07 + REPORT.md),
  `vignettes/articles/divas_case_study.Rmd` (precomputed article)
- PR #5: publication prep (v0.2.0) + this benchmark, on top of merged PR #4

## Caveats

- Concordance computed on scores (DIVAS RDS exposes no per-block
  factorization internals beyond loadings).
- cellprograms run without covariates for DIVAS parity; covariate handling
  was stress-tested separately (masked-refit analysis: classification flags
  78%/73% → 0% after artifact-gene masking).
- The severity-only baseline's high R² reflects wos being nearly determined
  by the 3-level severity band — it bounds both methods, not a flaw.
- Score concordance does not imply loading agreement: mutual-best pairs
  share sample-level axes (median top-100 loading Jaccard 0.16), so
  gene-level interpretations should be made within one method at a time.
