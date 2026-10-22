<div class="page-header">

# Summary of programs per cell type (n, K, R2 of the reconstruction, ELBO)

<div class="d-none name">

`program_summary.Rd`

</div>

</div>

<div class="ref-description section level2">

Summary of programs per cell type (n, K, R2 of the reconstruction, ELBO)

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
program_summary(fit)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
program_summary(fit)
```

</div>
