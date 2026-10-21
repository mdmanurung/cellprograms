# Program metadata table (one row per program)

Program metadata table (one row per program)

## Usage

``` r
program_metadata(fit)
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
program_metadata(fit)
#> Error in program_metadata(fit): could not find function "program_metadata"
```
