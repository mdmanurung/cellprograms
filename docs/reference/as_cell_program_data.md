<div class="page-header">

# Create a cell program data object

<div class="d-none name">

`as_cell_program_data.Rd`

</div>

</div>

<div class="ref-description section level2">

Create a cell program data object

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
as_cell_program_data(
  pseudobulk,
  sample_metadata = NULL,
  feature_metadata = NULL
)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

pseudobulk<a href="#arg-pseudobulk" class="anchor" aria-label="anchor"></a><br>
Named list of numeric matrices (observations x genes), one per cell type. Row names must be observation IDs.

sample_metadata<a href="#arg-sample-metadata" class="anchor" aria-label="anchor"></a><br>
A data.frame containing an \`observation_id\` column. Outcome variables may be present but are never used for factorization.

feature_metadata<a href="#arg-feature-metadata" class="anchor" aria-label="anchor"></a><br>
Optional data.frame keyed by gene ID.

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

An object of class \`cell_program_data\`.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
x
```

</div>
