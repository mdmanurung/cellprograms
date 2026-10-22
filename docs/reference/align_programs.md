<div class="page-header">

# Align all cell types' program subspaces to a reference cell type

<div class="d-none name">

`align_programs.Rd`

</div>

</div>

<div class="ref-description section level2">

Align all cell types' program subspaces to a reference cell type

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
align_programs(fit, reference = NULL, space = c("scores", "loadings"), ...)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_fit\` object.

reference<a href="#arg-reference" class="anchor" aria-label="anchor"></a><br>
Reference cell type (default: most programs).

space<a href="#arg-space" class="anchor" aria-label="anchor"></a><br>
Passed to \`principal_angles\` (default \`"scores"\`).

...<a href="#arg--" class="anchor" aria-label="anchor"></a><br>
Passed to \`principal_angles\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
al <- align_programs(fit, n_perm = 99, seed = 1)
attr(al, "reference")
al[[1]]$cosines
```

</div>
