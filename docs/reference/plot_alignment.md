<div class="page-header">

# Plot program alignment (principal angles)

<div class="d-none name">

`plot_alignment.Rd`

</div>

</div>

<div class="ref-description section level2">

Visualize \`principal_angles\` / \`align_programs\` results. For a single \`principal_angles\` object: how each program of each cell type loads on the shared principal directions (program-weight tiles), annotated with the permutation-based subspace test. For an \`align_programs\` result (named list of such objects): faceted per cell-type pair.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_alignment(x, q_alpha = 0.05)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

x<a href="#arg-x" class="anchor" aria-label="anchor"></a><br>
A \`principal_angles\` object, or the named list returned by \`align_programs\`.

q_alpha<a href="#arg-q-alpha" class="anchor" aria-label="anchor"></a><br>
Significance level for the per-pair subspace test (default 0.05).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A ggplot object.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
plot_alignment(principal_angles(fit, "B", "CD4", n_perm = 200))
```

</div>
