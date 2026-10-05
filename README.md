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

fit <- fit_celltype_programs(x, max_factors = 30)
fit <- canonicalize_programs(fit)

Z <- program_scores(fit)   # sample x cellular-program activity matrix
W <- program_loadings(fit) # gene x program loadings
```

## Status

v0.1.0 — core data model, flashier backend, canonicalization, and simulation
utilities. Stability assessment, integration backends (ICA, MOFA-FLEX,
EV-BIDIFAC), and plotting are on the roadmap (see `docs/final_report.md`).

Default configuration (`point_laplace` gene prior, per-gene variance) was
selected by a benchmark-first optimization loop against simulation
ground-truth recovery; see `docs/plan_review.md` and `docs/final_report.md`.

## License

MIT
