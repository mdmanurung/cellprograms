<div class="page-header">

# Plot program stability

<div class="d-none name">

`plot_stability.Rd`

</div>

</div>

<div class="ref-description section level2">

Dot plot of bootstrap stability per program: recovery frequency (fraction of bootstrap refits in which the program reappears) against the median score correlation between the original and matched bootstrap programs.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_stability(fit)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` object processed by \`assess_program_stability\` (uses its \`\$stability\` element).

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
fit <- assess_program_stability(fit, n_boot = 20, seed = 1)
plot_stability(fit)
```

</div>
