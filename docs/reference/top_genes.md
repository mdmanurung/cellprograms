# Top genes of a program

Top genes of a program

## Usage

``` r
top_genes(fit, program, n = 20L, by = c("abs_loading", "lfsr"))
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- program:

  Program ID, e.g. \`"Mono_3"\`.

- n:

  Number of genes.

- by:

  Rank by \`"abs_loading"\` or \`"lfsr"\` (flashier local false sign
  rate; needs the fit's gene-side lfsr).

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
top_genes(fit, "B_1", n = 5)
#>   gene   loading         lfsr
#> 2 G117  2.925181 2.920352e-40
#> 3  G41  2.824554 1.119050e-37
#> 4 G292  2.722820 3.720503e-35
#> 7  G64 -2.270005 0.000000e+00
#> 6 G130  2.252322 1.108095e-24
```
