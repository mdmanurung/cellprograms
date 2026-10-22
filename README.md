# cellprograms

[![R-CMD-check](https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdmanurung/cellprograms/actions/workflows/R-CMD-check.yaml)

Discover donor-level, cell-type-specific transcriptional programs from
pseudobulk expression using empirical-Bayes matrix factorization. Fit locally
within each cell type; investigate cross-type coordination afterwards rather
than assuming it.

**Development software, not a scientific release certificate.** The supported
primary route is ordinary local EBMF with explicit preprocessing, named axes,
effective estimator settings and failure accounting. Experimental integrations
and historical benchmark summaries are not established by ordinary-core tests.
See [support and migration](IMPLEMENTATION_SUPPORT.md) and the installed
`validation/CORE_SUPPORT.md` record.

## Installation

```r
remotes::install_github("mdmanurung/cellprograms")
```

Requires R >= 4.3. GitHub dependency revisions are pinned in `DESCRIPTION`.
For reproducibility, install a **full reviewed commit SHA**, not moving `main`,
and use the exact-source build/install receipt tooling in
[cellprograms-reproducibility](https://github.com/mdmanurung/cellprograms-reproducibility).
A package version string alone does not identify installed source.

## Ordinary local fit

```r
library(cellprograms)

# One observations x genes matrix per cell type, with unique named axes.
# metadata contains observation_id and any explicitly selected covariates.
x <- as_cell_program_data(pseudobulk, sample_metadata = metadata)
validate_cell_program_data(x)

fit <- fit_celltype_programs(
  x, max_factors = 10, feature_scale = "logcounts", seed = 1
)
fit <- canonicalize_programs(fit)

Z <- program_scores(fit)
W <- program_loadings(fit)
program_summary(fit)
fit$fits[[1]]$effective_config
```

Declare units truthfully; the software cannot infer that a matrix represents
counts, log counts, normalized expression or some other scale. Choice of rank,
preprocessing and covariates is analysis-specific, not a universal recommendation
established by the default arguments.

The normal activity / point-Laplace loading estimator inherits its native
flashier/ebnm behavior. Fallback attempts are warned and recorded; inspect
`effective_config` and `fallback`. Numeric `var_floor` is a **fixed known
residual standard deviation**, not a lower bound on an estimated variance.
Known-noise fits use a binding native path; they must not silently ignore it.

## Frozen preprocessing and projection

```r
spec <- freeze_preprocessing(
  x, features = "variable", n_variable_genes = 2000,
  feature_scale = "logcounts"
)
fit <- canonicalize_programs(fit_celltype_programs(
  x, preprocessing_spec = spec, max_factors = 10, seed = 1
))
projected <- project_programs_to_cells(
  fit, expression = query_matrix, cell_types = query_cell_types,
  feature_scale = "logcounts"
)
```

`query_matrix` has query observations in rows and genes in columns;
`query_cell_types` is named by query observation IDs. Supply donor/sample
metadata explicitly when the trained design requires it.

Frozen specifications verify their digest, feature identity, transformations
and units. Trained covariate transformations are replayed, not learned again
from the query composition. Query-time `feature_scale` describes **query** units;
it cannot relabel a reference with unknown units. Unknown-unit references,
legacy residualized fits without trained replay, or missing native manifests
require a truthful refit. Current round-trip tests do not certify every
historical package-version pair.

## Stability and honest denominators

```r
fit <- assess_program_stability(fit, n_boot = 20, seed = 1, progress = FALSE)
fit$stability$replicate_log
fit$stability$per_program
```

Inspect requested, successful and failed refits separately. The legacy
`recovery_freq` is `n_hit / n_boot`; `conditional_recovery_freq` uses successful
refits, including successful zero-rank fits. For a rate conditional on a
matchable cell type, use `n_hit / n_eligible` when positive. All-failure rates
are undefined (`NA`), not zero. Stability evaluates the full fitting procedure,
including its documented fallback ladder, not a guarantee of the same selected
prior on every resample.

## Coordination and sharing: interpretation limits

`find_coordinated_programs()` computes fixed-input permutation-max-adjusted
**node** p-values. Programs are not refitted under permutation. Validity depends
on the declared null, exchangeability, missingness and selection design; this
does not establish general learned-factor, edge-family or merged-CMP-set error
control. Ordinary row shuffling is not patient-cluster permutation for repeated
visits. Preserve uncallable nodes and reasons rather than interpreting them as
negative evidence.

`program_sharing()` and block memberships are **descriptive**. Activity in both
conditions does not prove equal effects, common structure, private biology or
calibrated inferential sharing. Prediction and invalid patient-pairing claims in
historical case studies remain withdrawn.

## Optional and historical surfaces

- `cebmf_torch` is a distinct experimental estimator, not numerical equivalence
  to ordinary flashier. Its supported controls and native environment must be
  validated independently; unsupported combinations fail early.
- SuperCell, SEACells/mcRigor, optional covariate integrations and metacell
  projection have separate dependency and aggregation-scale limitations.
  Ordinary checks do not certify them.
- Historical COVID, DIVAS, immune-aging and comparator material is retained for
  provenance. Maintained reconstruction and study status belong in the
  reproducibility repository. Neither a historical report nor an executed tiny
  fixture is new biological, calibration or independent-reproduction evidence.
- The [previous README](reports/historical/README_before_completion.md) is
  preserved as an explicitly superseded historical record, not current claims.

## Validation and documentation

```sh
Rscript -e 'testthat::test_local(".", stop_on_failure = TRUE)'
R CMD build .
# Use the actual tarball produced above, not an assumed release filename.
R CMD check --no-manual cellprograms_<version>.tar.gz
```

The installed test entrypoint executes ordinary native-fit and actual-refit
regressions as well as algebra, design replay, units, support and failure tests.
Expected warnings are asserted; unexpected or missing messages fail their tests.
Optional tests must identify their applicability instead of becoming empty
mandatory checks. CI defines Linux current/minimum-R, macOS and Windows lanes;
only actual run results establish coverage.

- [Quickstart](vignettes/quickstart.Rmd)
- [NEWS](NEWS.md)
- [Support and migration](IMPLEMENTATION_SUPPORT.md)
- [Reproducibility repository](https://github.com/mdmanurung/cellprograms-reproducibility)

## License and citation

MIT for the package's licensed code; preserve third-party notices. Dataset,
comparator and archived-result rights are separate. Use `citation("cellprograms")`
for package attribution and cite inherited methods separately.
