# Plot top-loaded genes per program

Bar plot of the largest-magnitude loadings for one or more programs,
colored by loading sign. Loadings are posterior means: small nonzero
values are shrinkage noise, so extract gene sets with an explicit
threshold or lfsr rather than exact zeros.

## Usage

``` r
plot_top_genes(fit, programs, n = 10)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\` object.

- programs:

  Program id(s) (e.g. \`"B_1"\` or \`c("B_1", "CD4_1")\`).

- n:

  Number of top genes per program (by \|loading\|).

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
plot_top_genes(fit, c("B_1", "CD4_1"), n = 8)
```
