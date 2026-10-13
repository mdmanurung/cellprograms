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
#>   program_id cell_type
#> 1        B_1         B
#> 2        B_2         B
#> 3      CD4_1       CD4
#> 4      CD8_1       CD8
#> 5       NK_1        NK
#> 6       NK_2        NK
#> 7     Mono_1      Mono
#> 8     Mono_2      Mono
```
