# cellprograms

An R package for discovering **cell-type-specific transcriptional programs** from
pseudobulk single-cell RNA-seq data using empirical-Bayes matrix factorization
(`flashier`), then optionally integrating those programs into higher-order
multicellular states.

**Discover, don't map.** Unlike reference-projection tools, `cellprograms`
identifies programs de novo within each cell type. The design principle is
*factorize locally first; integrate later only when scientifically useful* —
a B-cell program is discoverable because it explains B-cell expression,
regardless of whether anything analogous happens in other cell types.
Cross-cell-type coordination is investigated subsequently, never assumed.

## Installation

```r
remotes::install_github("mdmanurung/cellprograms")
```

Requires R >= 4.1. The `flashier` and `ebnm` backends are installed
automatically (`flashier` via GitHub).

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

## Status

v0.1.0 — core data model, flashier backend with hard (residualize) and soft
(fixed) covariate adjustment, canonicalization, a calibrated cross-cell-type
sharing test (principal angles in scores space, optional blocked-permutation
`strata`), subsampling stability, MOFAcellulaR export, and simulation
utilities. Integration backends such as ICA and EV-BIDIFAC were evaluated and
rejected; coordinated-program detection across cell types is in development
(see `ROADMAP.md`).

Default configuration (`point_laplace` gene prior, per-donor residual
variance, `max_factors = 10`) was selected by a benchmark-first optimization
loop against simulation ground-truth recovery and refined by a pseudobulk
tuning benchmark (below); see `docs/plan_review.md` and
`docs/final_report.md`.

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

## License

MIT
