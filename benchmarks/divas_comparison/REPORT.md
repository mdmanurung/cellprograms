# DIVAS COVID-19 case study: cellprograms reanalysis and comparison

Reanalysis of the DIVAS COVID-19 multi-omics case study
(https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy) with `cellprograms`,
and a systematic comparison of the two pipelines' outputs.

## Data

| Item | Value |
|---|---|
| Samples | 114 (85 patients: 29 with paired T1/T2, 56 with a single timepoint) |
| Timepoints | 57 T1 + 57 T2 samples |
| Severity | 40 mild / 40 moderate / 34 severe |
| Outcome | `wos` (WHO ordinal score) |
| Blocks | 4 gene blocks (CD4 T, CD8 T, CD14 monocytes, NK; 8,634 shared genes each), proteomics (481 proteins), metabolomics (763 metabolites) |

DIVAS was run by the original authors (`DIVASmain(nsim=500, colCent=FALSE,
rowCent=TRUE, seed=123)`) and its precomputed results were downloaded from
Zenodo (record 17430294). cellprograms was run at DIVAS parity
(`center=TRUE, scale=FALSE`, no covariates, `max_factors=30`,
`loading_prior="point_laplace"`, `var_type=1`, backfitting + nullcheck,
`seed=1`), which retained K=30 programs in every block.

## Pipeline

Scripts in `scripts/` (run in order):

1. `01_download_data.sh` - fetch the 6 block matrices, sample metadata, and
   the precomputed DIVAS RDS from Zenodo.
2. `02_fit_cellprograms.R` - fit cellprograms on all 6 blocks.
3. `03_cmp_sharing.R` - cross-block coordination (marginal and
   within-severity) and program sharing across the 4 gene blocks.
4. `04_extract_divas.R` - parse the DIVAS RDS into tidy tables
   (components, scores, loadings, block patterns).
5. `05_compare.R` - score concordance, taxonomy agreement, clinical
   associations, cross-validated `wos` prediction, unique findings.
6. `06_figures_extra.R` - concordance heatmap, breadth plot, effect-size
   agreement plot.
7. `07_sensitivity_pairing.R` - T1-only coordination and patient-blocked
   permutation null (pairing robustness).

## Results

### Program-level concordance

23 of 91 DIVAS components have a cellprograms counterpart at
|Spearman| > 0.5 (score correlation across the 114 samples); in the reverse
direction, 51 of 180 cellprograms scores have a DIVAS counterpart at
|Spearman| > 0.5. 42 component-program pairs are mutual best matches
(19 of them above 0.5). Every one of the 120 gene-block programs is
classified **shared** across cell types by the mashr analysis, and their
best-matching DIVAS components span all taxonomy classes: 96 partial,
18 individual, and all 6... (6 joint; `taxonomy_agreement.csv`). The single
DIVAS joint component has a strong best match (CD4_5, |rho| = 0.78) but is
not mutual-best.

Concordance increases with DIVAS pattern breadth: the median best-match
correlation rises from 0.38 (components private to 1 block) through 0.46-0.47
(2-3 blocks) and 0.56 (4 blocks) to 0.80 (5 blocks); the joint component
reaches 0.78 (`fig_breadth_concordance.png`). Cross-block structure - the
part of the biology DIVAS highlights - is recovered at least as well by
cellprograms' shared programs.

Score-level agreement does not extend to gene-level loadings: for the 40
mutual-best gene-block pairs, the top-100 loading Jaccard is low (median
0.16, IQR 0.14-0.22; `loading_overlap_top100.csv`) - the two factorizations
recover the same sample-level axes but assign genes to programs
differently.

Unique structure: 5 DIVAS components have no counterpart (all individual,
best |cor| < 0.3), and 13 cellprograms programs have no DIVAS counterpart -
cellprograms retains more block-specific factors (K=30 per block vs DIVAS's
91 total components shared across blocks).

### Clinical association

At BH < 0.05, 44/180 cellprograms scores associate with `wos` versus 8/91
DIVAS components (`clinical_associations.csv`). Cross-validated ridge
prediction of `wos` (`cv_wos_prediction.csv`):

| Predictors | n | R2 | RMSE |
|---|---|---|---|
| DIVAS components | 91 | 0.42 | 1.66 |
| cellprograms scores | 180 | 0.64 | 1.28 |
| severity-only baseline | 2 | 0.90 | 0.62 |

Both methods' scores are dominated by the severity axis (`wos` is nearly
determined by the 3-level severity band); the meaningful comparison is
between methods, where cellprograms scores capture ~50% more `wos`
variance than DIVAS components (0.64 vs 0.42).

For the 91 matched pairs, the agreement in clinical effect sizes is weak
overall (Spearman rho = 0.17), but pairs where either method is significant
at BH < 0.05 sit mostly above the diagonal: the cellprograms counterpart
carries as much or more `wos` signal than the DIVAS component
(`fig_effect_agreement.png`).

### Coordination and sharing

- Marginal coordination: 148 significant cross-block edges collapsed into
  6 coordinated program modules (CMPs).
- Within-severity coordination: 123 edges, 6 CMPs - the coordination
  structure is not an artifact of severity composition.
- Program sharing (mashr, 4 gene blocks): 120 programs in 31 groups
  (27 groups contain all 4 cell types); all 120 programs classify as
  shared (frac_top_shared >= 0.98).

The 120/120 shared result is more permissive than DIVAS's taxonomy
(44 individual / 46 partial / 1 joint). Two reasons: (i) the 4 immune cell
types in this cohort are highly correlated - DIVAS itself finds 46/91
components spanning >= 2 gene blocks - and the dominant axis (severity) is
present in every cell type; (ii) mash's lfsr-based classification is
liberal when conditions are strongly correlated. The two methods therefore
agree that cross-cell-type structure dominates this cohort, with
cellprograms attributing essentially all gene-block structure to shared
programs.

### Pairing robustness (T1-only + patient-blocked permutation)

Because 29 patients contribute two timepoints, the headline coordination
calls were checked two ways (`07_sensitivity_pairing.R`):

- **T1-only refit** (57 T1 samples, no repeated measures): 50/148
  full-cohort edges reproduce (34%), and both runs find 6 CMPs - the
  module structure is stable while the edge list shrinks with the halved
  sample size.
- **Patient-blocked permutation null** (all samples of a patient permuted
  together, 999 permutations): 121/148 edges (82%) remain significant at
  p < 0.05 (`sensitivity_patient_blocked_perm.csv`). The coordination
  structure is largely robust to the repeated-measures concern; the 18%
  of edges that do not survive were likely driven by patient-level
  composition.

## Caveats

- The DIVAS RDS provides component scores and loadings but not the exact
  per-block factorization; concordance is computed on scores.
- cellprograms was run without covariates for DIVAS parity; the
  publication thread's masked-refit analysis (separate) shows covariate
  handling affects program classification in other cohorts.
- The severity-only baseline's high R2 reflects the near-deterministic
  relationship between `wos` and the 3-level severity band, not a
  limitation of either method.
