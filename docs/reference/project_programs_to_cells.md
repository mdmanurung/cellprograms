<div class="page-header">

# Project fitted programs onto individual cells

<div class="d-none name">

`project_programs_to_cells.Rd`

</div>

</div>

<div class="ref-description section level2">

Estimates per-cell program activities for a \`cell_program_fit\`: each cell's expression vector is adjusted with the fit's recorded preprocessing (covariate residualization coefficients, then the stored centering and scale values), and activities are inferred against the frozen gene loadings with \`flashier.utils::project_factor_activity()\`. Raw activities are canonicalized with the same RMS-scale and sign multipliers as the donor scores, so projected activities are on the canonical program scale and comparable across cell types and cohorts.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
project_programs_to_cells(
  fit,
  expression,
  cell_types,
  feature_scale = NULL,
  donors = NULL,
  sample_metadata = NULL,
  min_coverage = 0.9,
  maxiter = 200L,
  query_var_type = 0L,
  dense_budget_bytes = 256 * 1024^2,
  keep_native = FALSE
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\` (from \`fit_celltype_programs()\` + \`canonicalize_programs()\`).

expression<a href="#arg-expression" class="anchor" aria-label="anchor"></a><br>
Dense numeric matrix with cells in rows and genes in columns (row and column names required). Sparse matrices are not supported in v1. Genes absent from a cell type's fit panel are handled by \`flashier.utils\` intersection matching against \`min_coverage\`.

cell_types<a href="#arg-cell-types" class="anchor" aria-label="anchor"></a><br>
Named character vector mapping every row of \`expression\` to its cell type; names must be exactly the row names of \`expression\`. Cell types not in the fit (and fit cell types with no cells here) yield \`NA\` activities.

feature_scale<a href="#arg-feature-scale" class="anchor" aria-label="anchor"></a><br>
Declared measurement units of \`expression\` (e.g. "logcounts"). Defaults to the fit's recorded \`feature_scale\`, which must not be "unknown"; a query whose declared units differ from the fit's is refused. This argument declares query units only: a reference with unknown or missing units must be refit with truthful units declared at fit time. It cannot be relabeled by a projection-time override.

donors<a href="#arg-donors" class="anchor" aria-label="anchor"></a><br>
Optional named character vector mapping every cell to its donor ID (\`sample_metadata\$observation_id\` values). Required whenever any projected cell type was fit with residualized covariates.

sample_metadata<a href="#arg-sample-metadata" class="anchor" aria-label="anchor"></a><br>
Optional data.frame with an \`observation_id\` column plus the covariate columns used at fit time. Required with \`donors\` when replaying residualization.

min_coverage<a href="#arg-min-coverage" class="anchor" aria-label="anchor"></a><br>
Minimum fraction of a cell type's fit genes that must be present in \`expression\` (passed to the projection oracle).

maxiter<a href="#arg-maxiter" class="anchor" aria-label="anchor"></a><br>
Maximum projection iterations per cell type.

query_var_type<a href="#arg-query-var-type" class="anchor" aria-label="anchor"></a><br>
Residual variance model for the query data: 0 = constant (default for single-cell profiles), 1 = per-cell, 2 = per-gene.

dense_budget_bytes<a href="#arg-dense-budget-bytes" class="anchor" aria-label="anchor"></a><br>
Memory budget for dense intermediates (passed to the projection oracle).

keep_native<a href="#arg-keep-native" class="anchor" aria-label="anchor"></a><br>
When TRUE, include the native flash fit object for each projected cell type in \`per_type\[\[ct\]\]\$native\`.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A list with components: \`activity\` (cells x programs matrix; canonicalized activities, \`NA\` for cells whose type contributed no programs), \`per_type\` (per cell type: \`activity\`, \`reconstruction\` metrics, and optionally \`native\`), \`coverage\` (per cell type: programs, cells, shared/reference genes, coverage), \`analysis_metadata\`, and \`provenance\`.

</div>

<div class="section level2">

## Details<a href="#details" class="anchor" aria-label="anchor"></a>

The fit must be canonicalized (\`canonicalize_programs()\`) and must carry a flashier.utils manifest (the default backend) with declared measurement units (\`feature_scale\`); legacy-backend and fixed-covariate fits are refused because their adjustment cannot be replayed. When the fit residualized covariates for a cell type, \`donors\` and \`sample_metadata\` are required so the adjustment can be rebuilt for the query cells; the rebuilt design must reproduce the fit's design columns exactly (fail closed on cohort changes that add or drop covariate levels).

Projection is joint over the supplied cells within each cell type: activity priors are re-estimated from the query set, so the same cell projected in two different query sets can receive different activities. Chunked projection would be an approximation and is not offered; oversize queries are refused pre-flight against \`dense_budget_bytes\` (decision record: flashier.utils \`reports/FU3B_CP6_CHUNKING_DECISION.md\`).

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable",
                             n_variable_genes = 150, feature_scale = "simulated")
fit <- canonicalize_programs(fit)
ct1 <- fit$cell_types[1]
# project the training pseudobulk back through the fit (donors as "cells")
cells <- setNames(rep(ct1, nrow(x$pseudobulk[[ct1]])), rownames(x$pseudobulk[[ct1]]))
proj <- project_programs_to_cells(fit, x$pseudobulk[[ct1]], cell_types = cells)
cor(proj$activity[, 1], fit$scores[[ct1]][, 1])   # high by construction
```

</div>
