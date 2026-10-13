# Create a cell program data object

Create a cell program data object

## Usage

``` r
as_cell_program_data(
  pseudobulk,
  sample_metadata = NULL,
  feature_metadata = NULL
)
```

## Arguments

- pseudobulk:

  Named list of numeric matrices (observations x genes), one per cell
  type. Row names must be observation IDs.

- sample_metadata:

  A data.frame containing an \`observation_id\` column. Outcome
  variables may be present but are never used for factorization.

- feature_metadata:

  Optional data.frame keyed by gene ID.

## Value

An object of class \`cell_program_data\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
x
#> cell_program_data
#> Cell types: 5 "B", "CD4", "CD8", "NK", and "Mono"
#> B: 24 observations x 300 genes
#> CD4: 24 observations x 300 genes
#> CD8: 24 observations x 300 genes
#> NK: 24 observations x 300 genes
#> Mono: 24 observations x 300 genes
#> Shared genes across all cell types: 300
```
