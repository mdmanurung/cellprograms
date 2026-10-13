# Program score matrix across all cell types

Program score matrix across all cell types

## Usage

``` r
program_scores(fit, format = c("wide", "long", "list"))
```

## Arguments

- fit:

  A canonicalized \`cell_program_fit\`.

- format:

  \`"wide"\` (observation x program, NA where a cell type was not
  observed; never imputed), \`"long"\`, or \`"list"\` (per-cell-type
  matrices).

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- fit_celltype_programs(x, max_factors = 3, features = "variable", n_variable_genes = 150)
fit <- canonicalize_programs(fit)
Z <- program_scores(fit)   # observations x programs, NA where unobserved
dim(Z)
#> [1] 24  8
head(program_scores(fit, "long"))
#>   observation_id cell_type program_id       score
#> 1           obs1         B        B_1  0.47142358
#> 2           obs2         B        B_1  0.06455346
#> 3           obs3         B        B_1 -0.83522990
#> 4           obs4         B        B_1  1.03868594
#> 5           obs5         B        B_1  0.91470905
#> 6           obs6         B        B_1  0.55071134
```
