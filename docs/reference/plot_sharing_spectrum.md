<div class="page-header">

# Plot the sharing spectrum

<div class="d-none name">

`plot_sharing_spectrum.Rd`

</div>

</div>

<div class="ref-description section level2">

Per-pair summary of \`sharing_spectrum\`: the smallest principal angle between each cell-type pair's program-activity subspaces, point size = number of shared subspace dimensions (cosine \> 0.5), color = permutation test outcome. Hollow points are underpowered pairs.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_sharing_spectrum(spectrum)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

spectrum<a href="#arg-spectrum" class="anchor" aria-label="anchor"></a><br>
Output of \`sharing_spectrum\`.

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
plot_sharing_spectrum(sharing_spectrum(fit, n_perm = 200))
```

</div>
