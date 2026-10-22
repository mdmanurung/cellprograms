<div class="page-header">

# Program score matrix across all cell types

<div class="d-none name">

`program_scores.Rd`

</div>

</div>

<div class="ref-description section level2">

Program score matrix across all cell types

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
program_scores(fit, format = c("wide", "long", "list"))
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

format<a href="#arg-format" class="anchor" aria-label="anchor"></a><br>
\`"wide"\` (observation x program, NA where a cell type was not observed; never imputed), \`"long"\`, or \`"list"\` (per-cell-type matrices).

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
Z <- program_scores(fit)   # observations x programs, NA where unobserved
dim(Z)
head(program_scores(fit, "long"))
```

</div>
