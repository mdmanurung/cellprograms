<div class="page-header">

# Validate a cell program data object

<div class="d-none name">

`validate_cell_program_data.Rd`

</div>

</div>

<div class="ref-description section level2">

Runs the full data-contract check and returns the object invisibly, emitting warnings for soft-failure conditions (small N, near-zero variance genes, extreme sparsity).

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
validate_cell_program_data(x)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

x<a href="#arg-x" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_data\` object.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
validate_cell_program_data(x)
```

</div>
