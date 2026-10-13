# cellprograms

[![R-CMD-check](https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml)

Discover cell-type-specific transcriptional programs from pseudobulk
single-cell RNA-seq, and integrate them into higher-order multicellular
states.

The cellprograms package fits empirical-Bayes matrix factorizations
(`flashier`) independently within each cell type, then — only when
scientifically useful — investigates how the recovered programs
coordinate across cell types, share their gene loadings, and associate
with context.

Major functionality of the cellprograms package:

- [`as_cell_program_data()`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)
  Build the validated pseudobulk data object (named list of observations
  x genes matrices, one per cell type)
- [`fit_celltype_programs()`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)
  Fit per-cell-type programs with the EBMF backend (hard residualization
  or soft fixed covariate adjustment)
- [`canonicalize_programs()`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)
  Align program sign, scale and labels across cell types
- [`program_scores()`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)
  Sample x program activity matrix across all cell types
- [`program_loadings()`](https://mdmanurung.github.io/cellprograms/reference/program_loadings.md)
  Gene x program loading matrix
- [`top_genes()`](https://mdmanurung.github.io/cellprograms/reference/top_genes.md)
  Top-loaded genes of a single program, with optional lfsr
- [`program_summary()`](https://mdmanurung.github.io/cellprograms/reference/program_summary.md)
  Per-cell-type summary: number of programs, R2, ELBO
- [`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
  Detect coordinated multicellular programs (CMPs) with family-wise
  error control by within-strata donor permutation
- [`cmp_scores()`](https://mdmanurung.github.io/cellprograms/reference/cmp_scores.md)
  Donor x CMP consensus scores
- [`annotate_cmp_tf()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_tf.md)
  TF footprint of every member program from CollecTRI regulons
  (competitive + self-contained tests)
- [`cmp_triage()`](https://mdmanurung.github.io/cellprograms/reference/cmp_triage.md)
  One-row-per-CMP triage: recurrent genes, nuisance counts, context
  association, replication, shared TFs
- [`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)
  Classify matched programs as shared, partial, or private (mashr on
  loading posteriors)
- [`artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/artifact_genes.md)
  /
  [`flag_artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md)
  Mask V(D)J and donor-genotype artifact genes
- [`simulate_pseudobulk()`](https://mdmanurung.github.io/cellprograms/reference/simulate_pseudobulk.md)
  Simulation with planted programs and ground truth

## Installation

\
`remotes``::`[`install_github`](https://remotes.r-lib.org/reference/install_github.html)`(``"mdmanurung/cellprograms"``)`

Requires R \>= 4.1. The `flashier` and `ebnm` backends are installed
automatically (`flashier` via GitHub).

## Quick start

\
[`library`](https://rdrr.io/r/base/library.html)`(``cellprograms``)`\
\
`# pseudobulk: named list of observations x genes matrices, one per cell type`\
`x`` ``<-`` `[`as_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)`(``pseudobulk``, sample_metadata ``=`` ``metadata``)`\
[`validate_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/validate_cell_program_data.md)`(``x``)`\
\
`fit`` ``<-`` `[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``x``)``  ``# max_factors = 10 default; override to revisit`\
`fit`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(``fit``)`\
\
`Z`` ``<-`` `[`program_scores`](https://mdmanurung.github.io/cellprograms/reference/program_scores.md)`(``fit``)``   ``# sample x cellular-program activity matrix`\
`W`` ``<-`` `[`program_loadings`](https://mdmanurung.github.io/cellprograms/reference/program_loadings.md)`(``fit``)`` ``# gene x program loadings`

## Motivation

**Discover, don’t map.** Unlike reference-projection tools,
`cellprograms` identifies programs de novo within each cell type. The
design principle is *factorize locally first; integrate later only when
scientifically useful* — a B-cell program is discoverable because it
explains B-cell expression, regardless of whether anything analogous
happens in other cell types. Cross-cell-type coordination is
investigated subsequently, never assumed.

## Resources

- Case studies: [real-data case
  study](https://mdmanurung.github.io/cellprograms/articles/real_data_case_study.html)
  (Stephenson + COMBAT COVID cohorts) and [DIVAS multi-omics
  reanalysis](https://mdmanurung.github.io/cellprograms/articles/divas_case_study.html)
- Pre-registered held-out benchmarks: `benchmarks/cmp_detection/` (CMP
  detection, 800 datasets per arm) and `benchmarks/divas_comparison/`
  (DIVAS comparison)
- Default-configuration tuning benchmark: see [Tuning
  benchmark](#tuning-benchmark-why-max_factors--10)

## Status

v0.2.0 — core data model, flashier backend with hard (residualize) and
soft (fixed) covariate adjustment, canonicalization, a calibrated
cross-cell-type sharing test (principal angles in scores space, optional
blocked-permutation `strata`), subsampling stability, MOFAcellulaR
export, and simulation utilities. Coordinated multicellular program
(CMP) detection is shipped and confirmed by two pre-registered held-out
benchmarks (`benchmarks/cmp_detection/`). Post-hoc shared-vs-specific
program assessment
([`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md),
mashr-based) and artifact-gene masking
([`artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/artifact_genes.md),
[`flag_artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md),
fit-time `exclude_genes`) are included. Integration backends such as ICA
and EV-BIDIFAC were evaluated and rejected (see `reports/ROADMAP.md`).

Default configuration (`point_laplace` gene prior, per-donor residual
variance, `max_factors = 10`) was selected by a benchmark-first
optimization loop against simulation ground-truth recovery and refined
by a pseudobulk tuning benchmark (below); see `reports/plan_review.md`
and `reports/final_report.md`.

## Tuning benchmark: why `max_factors = 10`?

The default configuration is grounded in a full-factorial tuning
benchmark on COMBAT COVID pseudobulk (770,833 cells, 134 donors, 16 cell
types), run with the same evaluators as the package
(`python/run_tuning_search.py` in this repository):

- **Search space**: 32 configs = loading prior (`point_laplace`,
  `point_normal`) x `max_factors` (10, 15, 30, 50) x `var_type` (0, 1) x
  backfit (TRUE, FALSE), on TMM + scran-HVG pseudobulk (edgeR TMM size
  factors, top-2000 HVGs per cell type).
- **Selection metric**: mean rank of disease-status kNN (corrected
  macro-F1) and WHO-ordinal severity kNN Spearman, with a reconstruction
  guardrail (held-out R^2 \>= 0.423, the MOFA2-comparator value).
- **Result**: `max_factors = 10` dominated every prior x variance x
  backfit combination. K \>= 30 overfits: with 16 cell types the
  concatenated representation reaches 400-600 factors and disease-signal
  recovery collapses to ~0. The winning config (`point_laplace`, K=10,
  `var_type = 1`, backfit) improved disease kNN 0.255 -\> 0.359 over the
  old K=30 default while passing the guardrail (held-out R^2 0.48-0.59
  across preprocessing arms).
- **Preprocessing parity**: the K-vs-biology relationship was stable
  across three arms (log-norm all genes, log-norm + HVG, TMM + HVG), so
  the default is not an artifact of one normalization.
- **`var_type`**: 1 (one residual-variance parameter per donor) is the
  right default for pseudobulk. `var_type = 2` (per-gene) crashes on
  log-normalized pseudobulk (NaN in the variance update) and is strictly
  worse where it runs. `c(1, 2)` (Kronecker rank-one) is supported but
  unbenchmarked at pseudobulk scale and substantially slower.
- **External validity**: on an independent cohort (Stephenson, 126
  donors, 43 cell types) the tuned config generalized (disease kNN 0.648
  vs 0.468 for the old default) and reduced retention of technical
  structure (pool-ID kNN ~0.46 -\> 0.21-0.24 on log-norm arms).

Benchmark artifacts (search tables, attribution, cross-dataset
validation, MOFA-FLEX comparison) and the summary report are available
on request; the runners live in `python/`.

## Coordinated multicellular programs

[`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
finds groups of programs from different cell types whose donor-level
scores co-vary (coordinated multicellular programs, CMPs), with
family-wise error control by within-strata donor permutation. Gene
overlap plays no role in the call;
[`annotate_cmp_lr()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_lr.md),
[`annotate_cmp_themes()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_themes.md),
[`cmp_gene_overlap()`](https://mdmanurung.github.io/cellprograms/reference/cmp_gene_overlap.md)
and
[`cmp_context_association()`](https://mdmanurung.github.io/cellprograms/reference/cmp_context_association.md)
annotate called CMPs, each with its own null.

\
`x`` ``<-`` `[`as_cell_program_data`](https://mdmanurung.github.io/cellprograms/reference/as_cell_program_data.md)`(``pseudobulk``, sample_metadata ``=`` ``metadata``)`\
`fit`` ``<-`` `[`canonicalize_programs`](https://mdmanurung.github.io/cellprograms/reference/canonicalize_programs.md)`(`[`fit_celltype_programs`](https://mdmanurung.github.io/cellprograms/reference/fit_celltype_programs.md)`(``x``, max_factors ``=`` ``30``, covariates ``=`` ``"site"``)``)`\
`res`` ``<-`` `[`find_coordinated_programs`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)`(``fit``, strata ``=`` ``"site"``, ambient ``=`` ``"flag"``, n_perm ``=`` ``999``)`\
[`cmp_scores`](https://mdmanurung.github.io/cellprograms/reference/cmp_scores.md)`(``fit``, ``res``)``                      ``# donors x CMPs consensus scores`

**Use `max_factors = 30` for fits intended for CMP detection.** A CMP
can only be called if each member program exists upstream. In the
pre-registered cmp-bench-v1.2 held-out benchmark (seeds 1051-1100; 8
scenarios x 2 simulation bases; 800 datasets per arm), fitting with
`max_factors = 30` raised pooled membership F1 from 0.571 to 0.705
(paired difference +0.127, 95% CI 0.113-0.142), entirely on the
semi-synthetic base built from Stephenson pseudobulk (0.671 -\> 0.934),
with every false-call guardrail \<= 0.06 (null FWER 0.04 / 0.00). The
global default stays at 10 because the tuning benchmark below optimizes
a different target (disease signal in the concatenated representation),
where K \>= 30 overfits.
[`find_coordinated_programs()`](https://mdmanurung.github.io/cellprograms/reference/find_coordinated_programs.md)
prints a one-time message when a fit with `max_factors < 30` has a cell
type at its cap. On real cohorts (Stephenson, COMBAT), cap 30 roughly
doubled the number of CMPs; many added CMPs carried immunoglobulin/TCR V
genes or genes with known donor-genotype variation (e.g. GSTM1,
HLA-DQ/DR haplotype genes). These may co-vary across cell types for
genetic rather than biological-coordination reasons; this is unconfirmed
without genotype data. Use `ambient = "flag"` and inspect each CMP’s top
genes,
[`cmp_context_association()`](https://mdmanurung.github.io/cellprograms/reference/cmp_context_association.md)
and replication before interpreting it. Benchmark scripts,
pre-registrations and result tables are in `benchmarks/cmp_detection/`.

### Interpreting CMPs: TF footprints and triage

\
`net`` ``<-`` `[`prepare_collectri`](https://mdmanurung.github.io/cellprograms/reference/prepare_collectri.md)`(``omnipath_collectri_table``)``   ``# OmniPath interactions?datasets=collectri`\
`tf``  ``<-`` `[`annotate_cmp_tf`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_tf.md)`(``fit``, ``cmp``, ``net``)``                ``# TF footprint of every member program`\
`tri`` ``<-`` `[`cmp_triage`](https://mdmanurung.github.io/cellprograms/reference/cmp_triage.md)`(``fit``, ``cmp``, context ``=`` ``"severity"``, covariates ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"institute"``, ``"sex"``)``,`\
`                  replicate_by ``=`` ``"institute"``, tf ``=`` ``tf``)`

[`annotate_cmp_tf()`](https://mdmanurung.github.io/cellprograms/reference/annotate_cmp_tf.md)
fits a univariate linear model of each member’s CMP-oriented loadings on
each signed CollecTRI regulon. A TF is called in a member only if it
passes both a competitive test (BH-adjusted ULM, targets vs other genes)
and a self-contained sign-flip test (targets coherently signed), which
stops regulons of unloaded genes from looking “depleted” in skewed
programs. TFs called with one sign in \>= 2 members and \>= half of the
members are `shared` (a candidate common upstream driver); `specific`
TFs show how one cell type reads out the signal. Footprints describe
which regulons a program resembles; they are not evidence that the TF
causes the coordination, and they can be circular (HLA-class-II CMPs
light up CIITA/RFX5 because the HLA genes *are* their targets).

[`cmp_triage()`](https://mdmanurung.github.io/cellprograms/reference/cmp_triage.md)
returns one row per CMP with its recurrent top genes, counts of nuisance
genes
([`cmp_nuisance_genes()`](https://mdmanurung.github.io/cellprograms/reference/cmp_nuisance_genes.md):
V(D)J, candidate genotype markers, sex chromosome, immediate-early
genes, platelet/erythrocyte markers), the ambient fraction, context
association, within-level sign replication and shared TFs, plus a
heuristic `label` for prioritization.

## Program sharing and artifact masking

[`program_sharing()`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)
assesses whether programs matched across cell types share their gene
loadings (mashr on the stored loading posteriors), classifying each
matched group as shared, partial, or private — correlated loadings alone
do not imply shared effects.
[`flag_artifact_genes()`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md)
marks programs whose top genes are V(D)J segments or known
donor-genotype markers;
`fit_celltype_programs(..., exclude_genes = artifact_genes(genes))`
removes them before variable-gene selection.

\
`res`` ``<-`` `[`program_sharing`](https://mdmanurung.github.io/cellprograms/reference/program_sharing.md)`(``fit``, n_perm ``=`` ``199``, z_thresh ``=`` ``1.5``)`\
`flags`` ``<-`` `[`flag_artifact_genes`](https://mdmanurung.github.io/cellprograms/reference/flag_artifact_genes.md)`(``fit``)`

On the Stephenson and COMBAT COVID cohorts (K = 30), the broadest
cross-cell-type “shared” program was the immunoglobulin/V(D)J artifact
program (ambient/clonal leakage), while interferon and cell-cycle
programs replicated as genuinely shared across both cohorts — mask
artifact genes before interpreting sharing.

## License

MIT
