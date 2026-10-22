<div class="page-header">

# Ambient-likeness of programs

<div class="d-none name">

`ambient_program_scores.Rd`

</div>

</div>

<div class="ref-description section level2">

For each program, the largest absolute correlation between its gene loadings and another cell type's expression specificity (\`profile_other - profile_own\`, log scale). High values indicate a program that tracks contamination by another cell type's transcripts (ambient RNA, doublets) rather than biology intrinsic to the cell type. Genuine biology can also score high (e.g. a monocyte-like program in DCs), so treat it as a flag.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
ambient_program_scores(fit, profiles = NULL, min_genes = 50L)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

fit<a href="#arg-fit" class="anchor" aria-label="anchor"></a><br>
A canonicalized \`cell_program_fit\`.

profiles<a href="#arg-profiles" class="anchor" aria-label="anchor"></a><br>
Optional named list of per-cell-type mean log-expression vectors named by gene; defaults to column means of \`fit\$pseudobulk\`.

min_genes<a href="#arg-min-genes" class="anchor" aria-label="anchor"></a><br>
Minimum shared genes to compare two cell types (default 50).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

Data frame with \`cell_type\`, \`program_id\`, \`ambient_score\` and \`source\` (the cell type giving the maximum).

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
ambient_program_scores(fit)
```

</div>
