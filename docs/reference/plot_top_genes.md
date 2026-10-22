<div class="page-header">

# Plot top-loaded genes per program

<div class="d-none name">

`plot_top_genes.Rd`

</div>

</div>

<div class="ref-description section level2">

Bar plot of the largest-magnitude loadings for one or more programs, colored by loading sign. Loadings are posterior means: small nonzero values are shrinkage noise, so extract gene sets with an explicit threshold or lfsr rather than exact zeros.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_top_genes(fit, programs, n = 10)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\` object.

programs<a href="#arg-programs" class="anchor" aria-label="anchor"></a><br>
Program id(s) (e.g. \`"B_1"\` or \`c("B_1", "CD4_1")\`).

n<a href="#arg-n" class="anchor" aria-label="anchor"></a><br>
Number of top genes per program (by \|loading\|).

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
plot_top_genes(fit, c("B_1", "CD4_1"), n = 8)
```

</div>
