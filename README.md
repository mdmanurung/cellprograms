# cellprograms

[![R-CMD-check](https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml)

Discover cell-type-specific transcriptional programs from pseudobulk single-cell RNA-seq, and integrate them into higher-order multicellular states.

The cellprograms package fits empirical-Bayes matrix factorizations
(`flashier`) independently within each cell type, then — only when
scientifically useful — investigates how the recovered programs coordinate
across cell types, share their gene loadings, and associate with context.

Major functionality of the cellprograms package:

- `as_cell_program_data()` Build the validated pseudobulk data object (named list of observations x genes matrices, one per cell type)
- `fit_celltype_programs()` Fit per-cell-type programs with the EBMF backend (hard residualization or soft fixed covariate adjustment; `flashier` or `cebmf_torch` backend)
- `freeze_preprocessing()` / `project_programs_to_cells()` Freeze a fit's preprocessing spec (MR-0) and project new samples onto its programs through a unit-checked oracle
- `assess_program_stability()` Subsampling stability of recovered programs, with explicit accounting of failed refits
- `canonicalize_programs()` Align program sign, scale and labels across cell types
- `program_scores()` Sample x program activity matrix across all cell types
- `program_loadings()` Gene x program loading matrix
- `top_genes()` Top-loaded genes of a single program, with optional lfsr
- `program_summary()` Per-cell-type summary: number of programs, R2, ELBO
- `find_coordinated_programs()` Detect coordinated multicellular programs (CMPs) with family-wise error control by within-strata donor permutation
- `cmp_scores()` Donor x CMP consensus scores
- `annotate_cmp_tf()` TF footprint of every member program from CollecTRI regulons (competitive + self-contained tests)
- `cmp_triage()` One-row-per-CMP triage: recurrent genes, nuisance counts, context association, replication, shared TFs
- `program_sharing()` Classify matched programs as shared, partial, or private (mashr on loading posteriors)
- `artifact_genes()` / `flag_artifact_genes()` Mask V(D)J and donor-genotype artifact genes
- `simulate_pseudobulk()` Simulation with planted programs and ground truth

## Installation

```r
remotes::install_github("mdmanurung/cellprograms")
```

Requires R >= 4.3. The `flashier` and `ebnm` backends are installed
automatically (`flashier` via GitHub). The alternative `cebmf_torch` backend
(optionally selected with `fit_celltype_programs(backend = "cebmf_torch")`)
runs in a pinned Python environment (torch 2.7.1 + `cebmf_torch`); set
`options(cellprograms.cebmf_python_env = ...)` to point at an existing
environment, otherwise a basilisk environment is provisioned.

## Quick start

```r
library(cellprograms)

# pseudobulk: named list of observations x genes matrices, one per cell type
x <- as_cell_program_data(pseudobulk, sample_metadata = metadata)
validate_cell_program_data(x)

fit <- fit_celltype_programs(x)  # max_factors = 10 default; override to revisit
fit <- canonicalize_programs(fit)

Z <- program_scores(fit)   # sample x cellular-program activity matrix
W <- program_loadings(fit) # gene x program loadings
```

## Frozen preprocessing and out-of-sample projection

Fits record their full preprocessing specification, and
`project_programs_to_cells()` refuses queries whose declared measurement
units differ from the source fit's — a projection oracle, not a silent
re-derivation:

```r
spec <- freeze_preprocessing(x, features = "variable", n_variable_genes = 2000)
fit  <- fit_celltype_programs(x, preprocessing_spec = spec,
                              feature_scale = "logcounts")
Z_new <- project_programs_to_cells(fit, new_pseudobulk,
                                   feature_scale = "logcounts")
```

A refit supplied with a frozen spec cannot silently drift from the frozen
panel: any explicit conflict with the spec is an error. Stability
(`assess_program_stability()`) reports how many subsample refits failed
(`n_failed`/`n_success`) so that recovery frequencies are computed over
successful refits, never silently diluted.

## Motivation

**Discover, don't map.** Unlike reference-projection tools, `cellprograms`
identifies programs de novo within each cell type. The design principle is
*factorize locally first; integrate later only when scientifically useful* —
a B-cell program is discoverable because it explains B-cell expression,
regardless of whether anything analogous happens in other cell types.
Cross-cell-type coordination is investigated subsequently, never assumed.

## Resources

- Getting started: [quickstart](https://mdmanurung.github.io/cellprograms/articles/quickstart.html)
- Methods: [program mechanisms and the data model](https://mdmanurung.github.io/cellprograms/articles/model_mechanisms.html)
- Case studies: [real-data case study](https://mdmanurung.github.io/cellprograms/articles/real_data_case_study.html) (Stephenson + COMBAT COVID cohorts) and [DIVAS multi-omics reanalysis](https://mdmanurung.github.io/cellprograms/articles/divas_case_study.html)
- Application: [immune-aging cross-cohort analysis](https://mdmanurung.github.io/cellprograms/articles/immune_aging_cross_cohort.html) (8 PBMC cohorts, 2,609 donors; conservation panel + composition axis; scripts and tables under `analyses/immune_aging/`)
- Pre-registered held-out benchmarks: `benchmarks/cmp_detection/` (CMP detection, 800 datasets per arm) and `benchmarks/divas_comparison/` (DIVAS comparison)
- Default-configuration tuning benchmark: see [Tuning benchmark](#tuning-benchmark-why-max_factors--10)

## Status

v0.2.0 — core data model, flashier backend with hard (residualize) and soft
(fixed) covariate adjustment, canonicalization, a calibrated cross-cell-type
sharing test (principal angles in scores space, optional blocked-permutation
`strata`), subsampling stability, MOFAcellulaR export, and simulation
utilities. Coordinated multicellular program (CMP) detection is shipped and
confirmed by two pre-registered held-out benchmarks
(`benchmarks/cmp_detection/`). Post-hoc shared-vs-specific program assessment
(`program_sharing()`, mashr-based) and artifact-gene masking
(`artifact_genes()`, `flag_artifact_genes()`, fit-time `exclude_genes`) are
included. Integration backends such as ICA and EV-BIDIFAC were evaluated and
rejected (see `reports/ROADMAP.md`).

Recent additions: an alternative `cebmf_torch` torch backend (same
preprocessing and output contract; fixed residual SDs bind there, while the
flashier path currently passes them through unenforced — see the
`var_floor` docs), frozen preprocessing specifications with a unit-checked
projection oracle (`freeze_preprocessing()`,
`project_programs_to_cells()`), and a hardening pass across the analysis
contract: stability refits account for failed subsamples explicitly, CMP
member orientation propagates along the linking edges that formed the CMP,
and `program_sharing()` validates `matches=` references and warns on forced
pairs ("shared" means active-in-both, not proven shared structure — check
`p_all_shared`).

Default configuration (`point_laplace` gene prior, per-donor residual
variance, `max_factors = 10`) was selected by a benchmark-first optimization
loop against simulation ground-truth recovery and refined by a pseudobulk
tuning benchmark (below); see `reports/plan_review.md` and
`reports/final_report.md`.

## Tuning benchmark: why `max_factors = 10`?

The default configuration is grounded in a full-factorial tuning benchmark on
COMBAT COVID pseudobulk (770,833 cells, 134 donors, 16 cell types), run with
the same evaluators as the package (`python/run_tuning_search.py` in this
repository):

- **Search space**: 32 configs = loading prior (`point_laplace`, `point_normal`)
  x `max_factors` (10, 15, 30, 50) x `var_type` (0, 1) x backfit (TRUE, FALSE),
  on TMM + scran-HVG pseudobulk (edgeR TMM size factors, top-2000 HVGs per
  cell type).
- **Selection metric**: mean rank of disease-status kNN (corrected macro-F1)
  and WHO-ordinal severity kNN Spearman, with a reconstruction guardrail
  (held-out R^2 >= 0.423, the MOFA2-comparator value).
- **Result**: `max_factors = 10` dominated every prior x variance x backfit
  combination. K >= 30 overfits: with 16 cell types the concatenated
  representation reaches 400-600 factors and disease-signal recovery collapses
  to ~0. The winning config (`point_laplace`, K=10, `var_type = 1`, backfit)
  improved disease kNN 0.255 -> 0.359 over the old K=30 default while passing
  the guardrail (held-out R^2 0.48-0.59 across preprocessing arms).
- **Preprocessing parity**: the K-vs-biology relationship was stable across
  three arms (log-norm all genes, log-norm + HVG, TMM + HVG), so the default
  is not an artifact of one normalization.
- **`var_type`**: 1 (one residual-variance parameter per donor) is the right
  default for pseudobulk. `var_type = 2` (per-gene) crashes on log-normalized
  pseudobulk (NaN in the variance update) and is strictly worse where it runs.
  `c(1, 2)` (Kronecker rank-one) is supported but unbenchmarked at pseudobulk
  scale and substantially slower.
- **External validity**: on an independent cohort (Stephenson, 126 donors,
  43 cell types) the tuned config generalized (disease kNN 0.648 vs 0.468 for
  the old default) and reduced retention of technical structure (pool-ID kNN
  ~0.46 -> 0.21-0.24 on log-norm arms).

Benchmark artifacts (search tables, attribution, cross-dataset validation,
MOFA-FLEX comparison) and the summary report are available on request; the
runners live in `python/`.

## Coordinated multicellular programs

`find_coordinated_programs()` finds groups of programs from different cell
types whose donor-level scores co-vary (coordinated multicellular programs,
CMPs), with family-wise error control by within-strata donor permutation.
Gene overlap plays no role in the call; `annotate_cmp_lr()`,
`annotate_cmp_themes()`, `cmp_gene_overlap()` and `cmp_context_association()`
annotate called CMPs, each with its own null.

```r
x <- as_cell_program_data(pseudobulk, sample_metadata = metadata)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 30, covariates = "site"))
res <- find_coordinated_programs(fit, strata = "site", ambient = "flag", n_perm = 999)
cmp_scores(fit, res)                      # donors x CMPs consensus scores
```

**Use `max_factors = 30` for fits intended for CMP detection.** A CMP can only
be called if each member program exists upstream. In the pre-registered
cmp-bench-v1.2 held-out benchmark (seeds 1051-1100; 8 scenarios x 2
simulation bases; 800 datasets per arm), fitting with `max_factors = 30`
raised pooled membership F1 from 0.571 to 0.705 (paired difference +0.127,
95% CI 0.113-0.142), entirely on the semi-synthetic base built from
Stephenson pseudobulk (0.671 -> 0.934), with every false-call guardrail
<= 0.06 (null FWER 0.04 / 0.00). The global default stays at 10 because the
tuning benchmark below optimizes a different target (disease signal in the
concatenated representation), where K >= 30 overfits.
`find_coordinated_programs()` prints a one-time message when a fit with
`max_factors < 30` has a cell type at its cap. On real cohorts (Stephenson,
COMBAT), cap 30 roughly doubled the number of CMPs; many added CMPs carried
immunoglobulin/TCR V genes or genes with known donor-genotype variation (e.g.
GSTM1, HLA-DQ/DR haplotype genes). These may co-vary across cell types for
genetic rather than biological-coordination reasons; this is unconfirmed
without genotype data. Use `ambient = "flag"` and inspect
each CMP's top genes, `cmp_context_association()` and replication before
interpreting it. Benchmark scripts,
pre-registrations and result tables are in `benchmarks/cmp_detection/`.

### Interpreting CMPs: TF footprints and triage

```r
net <- prepare_collectri(omnipath_collectri_table)   # OmniPath interactions?datasets=collectri
tf  <- annotate_cmp_tf(fit, cmp, net)                # TF footprint of every member program
tri <- cmp_triage(fit, cmp, context = "severity", covariates = c("institute", "sex"),
                  replicate_by = "institute", tf = tf)
```

`annotate_cmp_tf()` fits a univariate linear model of each member's
CMP-oriented loadings on each signed CollecTRI regulon. A TF is called in a
member only if it passes both a competitive test (BH-adjusted ULM, targets vs
other genes) and a self-contained sign-flip test (targets coherently signed),
which stops regulons of unloaded genes from looking "depleted" in skewed
programs. TFs called with one sign in >= 2 members and >= half of the members
are `shared` (a candidate common upstream driver); `specific` TFs show how one
cell type reads out the signal. Footprints describe which regulons a program
resembles; they are not evidence that the TF causes the coordination, and they
can be circular (HLA-class-II CMPs light up CIITA/RFX5 because the HLA genes
*are* their targets).

`cmp_triage()` returns one row per CMP with its recurrent top genes, counts of
nuisance genes (`cmp_nuisance_genes()`: V(D)J, candidate genotype markers, sex
chromosome, immediate-early genes, platelet/erythrocyte markers), the ambient
fraction, context association, within-level sign replication and shared TFs,
plus a heuristic `label` for prioritization.

## Program sharing and artifact masking

`program_sharing()` assesses whether programs matched across cell types share
their gene loadings (mashr on the stored loading posteriors), classifying each
matched group as shared, partial, or private — correlated loadings alone do
not imply shared effects. `flag_artifact_genes()` marks programs whose top
genes are V(D)J segments or known donor-genotype markers;
`fit_celltype_programs(..., exclude_genes = artifact_genes(genes))` removes
them before variable-gene selection.

```r
res <- program_sharing(fit, n_perm = 199, z_thresh = 1.5)
flags <- flag_artifact_genes(fit)
```

On the Stephenson and COMBAT COVID cohorts (K = 30), the broadest cross-cell-type
"shared" program was the immunoglobulin/V(D)J artifact program (ambient/clonal
leakage), while interferon and cell-cycle programs replicated as genuinely
shared across both cohorts — mask artifact genes before interpreting sharing.

## License

MIT
