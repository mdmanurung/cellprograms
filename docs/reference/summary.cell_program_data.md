<div class="page-header">

# Summarize a cell program data object

<div class="d-none name">

`summary.cell_program_data.Rd`

</div>

</div>

<div class="ref-description section level2">

Summarize a cell program data object

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
# S3 method for class 'cell_program_data'
summary(object, ...)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

object<a href="#arg-object" class="anchor" aria-label="anchor"></a><br>
A \`cell_program_data\` object.

...<a href="#arg--" class="anchor" aria-label="anchor"></a><br>
Ignored.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
summary(x)
```

</div>
