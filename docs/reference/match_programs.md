<div class="page-header">

# Match programs across two fits (bootstrap refits or cohorts)

<div class="d-none name">

`match_programs.Rd`

</div>

</div>

<div class="ref-description section level2">

Match programs across two fits (bootstrap refits or cohorts)

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
match_programs(
  fit_a,
  fit_b,
  cell_type,
  method = c("subspace", "correlation"),
  space = c("loadings", "scores"),
  sim_threshold = 0.5
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit_a, fit_b<a href="#arg-fit-a-fit-b" class="anchor" aria-label="anchor"></a><br>
\`cell_program_fit\` objects.

cell_type<a href="#arg-cell-type" class="anchor" aria-label="anchor"></a><br>
Cell type to match (must exist in both fits).

method<a href="#arg-method" class="anchor" aria-label="anchor"></a><br>
\`"subspace"\` (default) or \`"correlation"\`. Factor-level \`matches\` are the same greedy loading-cosine matches in both; \`"subspace"\` additionally returns principal \`angles\` and a \`split_map\` (heuristic thresholds: cosine \> 0.5, \|weight\| \> 0.3; weights carry a singular-value scale and are not comparable across cell types).

space<a href="#arg-space" class="anchor" aria-label="anchor"></a><br>
Geometry for the subspace \`angles\`/\`split_map\`, used only when \`method = "subspace"\`: \`"loadings"\` (default; gene identity) or \`"scores"\` (activity; requires \>= 2 shared observations between the fits). The factor-level \`matches\` are always loading-cosine based, regardless of \`space\`.

sim_threshold<a href="#arg-sim-threshold" class="anchor" aria-label="anchor"></a><br>
Minimum \|cosine\| to report a factor-level match.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
# a second fit of the same data (different seed) stands in for a bootstrap refit
fit2 <- fit_celltype_programs(x, max_factors = 3, features = "variable",
                              n_variable_genes = 150, seed = 2)
fit2 <- canonicalize_programs(fit2)
match_programs(fit, fit2, "B")$matches
```

</div>
