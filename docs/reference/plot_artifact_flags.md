<div class="page-header">

# Plot artifact-gene flags

<div class="d-none name">

`plot_artifact_flags.Rd`

</div>

</div>

<div class="ref-description section level2">

Bar plot of the fraction of each program's top genes that are V(D)J or donor-genotype-marker genes, with the flag threshold marked.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_artifact_flags(flags, threshold = 0.2)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

flags<a href="#arg-flags" class="anchor" aria-label="anchor"></a><br>
Output of \`flag_artifact_genes\`.

threshold<a href="#arg-threshold" class="anchor" aria-label="anchor"></a><br>
Flag threshold (drawn as a dashed line; default 0.2).

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
plot_artifact_flags(flag_artifact_genes(fit))
```

</div>
