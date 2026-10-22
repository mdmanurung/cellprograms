<div class="page-header">

# Program loading matrices

<div class="d-none name">

`program_loadings.Rd`

</div>

</div>

<div class="ref-description section level2">

Program loading matrices

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
program_loadings(fit, format = c("list", "long"))
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

format<a href="#arg-format" class="anchor" aria-label="anchor"></a><br>
\`"list"\` (gene x program per cell type) or \`"long"\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
W <- program_loadings(fit)
sapply(W, dim)
```

</div>
