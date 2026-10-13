# Align all cell types' program subspaces to a reference cell type

Align all cell types' program subspaces to a reference cell type

## Usage

``` r
align_programs(fit, reference = NULL, space = c("scores", "loadings"), ...)
```

## Arguments

- fit:

  A \`cell_program_fit\` object.

- reference:

  Reference cell type (default: most programs).

- space:

  Passed to \`principal_angles\` (default \`"scores"\`).

- ...:

  Passed to \`principal_angles\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
al <- align_programs(fit, n_perm = 99, seed = 1)
attr(al, "reference")
#> [1] "B"
al[[1]]$cosines
#> [1] 0.9803336
```
