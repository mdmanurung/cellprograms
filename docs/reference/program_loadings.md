# Program loading matrices

Program loading matrices

## Usage

``` r
program_loadings(fit, format = c("list", "long"))
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- format:

  \`"list"\` (gene x program per cell type) or \`"long"\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
W <- program_loadings(fit)
sapply(W, dim)
#>        B CD4 CD8  NK Mono
#> [1,] 150 150 150 150  150
#> [2,]   2   1   1   2    2
```
