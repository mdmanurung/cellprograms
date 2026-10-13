# Summary of programs per cell type (n, K, R2 of the reconstruction, ELBO)

Summary of programs per cell type (n, K, R2 of the reconstruction, ELBO)

## Usage

``` r
program_summary(fit)
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
program_summary(fit)
#>   cell_type n_obs n_genes n_programs        R2      elbo
#> 1         B    24     150          2 0.4107044 -5719.888
#> 2       CD4    24     150          1 0.3463505 -5544.727
#> 3       CD8    24     150          1 0.2399331 -5506.694
#> 4        NK    24     150          2 0.3943826 -5716.026
#> 5      Mono    24     150          2 0.4113416 -5653.769
```
