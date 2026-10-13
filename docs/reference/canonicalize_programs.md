# Canonicalize program scale and sign

RMS-scales sample scores to unit root-mean-square with the compensating
factor applied to gene loadings (reconstruction preserved exactly), and
orients each factor so its largest-\|loading\| gene is positive.

## Usage

``` r
canonicalize_programs(fit)
```

## Arguments

- fit:

  A \`cell_program_fit\` object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
head(fit$scores$B[, 1])     # canonicalized donor scores of program B_1
#>        obs1        obs2        obs3        obs4        obs5        obs6 
#>  0.47142358  0.06455346 -0.83522990  1.03868594  0.91470905  0.55071134 
head(fit$loadings$B[, 1])   # canonicalized gene loadings of B_1
#>         G160         G117          G41         G292         G172         G130 
#> -0.010419810  2.925181026  2.824553881  2.722820113 -0.001974884  2.252321560 
```
