# Plot variance explained per program

Tile plot of the fraction of variance each program explains in its own
cell type's pseudobulk matrix (from the underlying flashier fit).

## Usage

``` r
plot_variance_explained(fit)
```

## Arguments

- fit:

  A \`cell_program_fit\` object.

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
plot_variance_explained(fit)
```
