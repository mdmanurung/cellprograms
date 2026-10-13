# Plot program stability

Dot plot of bootstrap stability per program: recovery frequency
(fraction of bootstrap refits in which the program reappears) against
the median score correlation between the original and matched bootstrap
programs.

## Usage

``` r
plot_stability(fit)
```

## Arguments

- fit:

  A \`cell_program_fit\` object processed by
  \`assess_program_stability\` (uses its \`\$stability\` element).

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
fit <- assess_program_stability(fit, n_boot = 20, seed = 1)
#> subsample 10/20
#> subsample 20/20
plot_stability(fit)
```
