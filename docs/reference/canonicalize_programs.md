<div class="page-header">

# Canonicalize program scale and sign

<div class="d-none name">

`canonicalize_programs.Rd`

</div>

</div>

<div class="ref-description section level2">

RMS-scales sample scores to unit root-mean-square with the compensating factor applied to gene loadings (reconstruction preserved exactly), and orients each factor so its largest-\|loading\| gene is positive.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
canonicalize_programs(fit)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` object.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
head(fit$scores$B[, 1])     # canonicalized donor scores of program B_1
head(fit$loadings$B[, 1])   # canonicalized gene loadings of B_1
```

</div>
