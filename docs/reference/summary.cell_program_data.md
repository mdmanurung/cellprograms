# Summarize a cell program data object

Summarize a cell program data object

## Usage

``` r
# S3 method for class 'cell_program_data'
summary(object, ...)
```

## Arguments

- object:

  A \`cell_program_data\` object.

- ...:

  Ignored.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
summary(x)
#> cell_program_data
#> Cell types: 5 "B", "CD4", "CD8", "NK", and "Mono"
#> B: 24 observations x 300 genes
#> CD4: 24 observations x 300 genes
#> CD8: 24 observations x 300 genes
#> NK: 24 observations x 300 genes
#> Mono: 24 observations x 300 genes
#> Shared genes across all cell types: 300
```
