<div class="page-header">

# Top genes of a program

<div class="d-none name">

`top_genes.Rd`

</div>

</div>

<div class="ref-description section level2">

Top genes of a program

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
top_genes(fit, program, n = 20L, by = c("abs_loading", "lfsr"))
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

program<a href="#arg-program" class="anchor" aria-label="anchor"></a><br>
Program ID, e.g. \`"Mono_3"\`.

n<a href="#arg-n" class="anchor" aria-label="anchor"></a><br>
Number of genes.

by<a href="#arg-by" class="anchor" aria-label="anchor"></a><br>
Rank by \`"abs_loading"\` or \`"lfsr"\` (flashier local false sign rate; needs the fit's gene-side lfsr).

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
top_genes(fit, "B_1", n = 5)
```

</div>
