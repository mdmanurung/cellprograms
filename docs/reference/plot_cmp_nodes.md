# Plot CMP node statistics

Bar plot of each program's node statistic A_v (pooled coordination
evidence across partner cell types), colored by CMP membership, with
permutation p-values marked.

## Usage

``` r
plot_cmp_nodes(cmp)
```

## Arguments

- cmp:

  Output of \`find_coordinated_programs\`.

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
cmp <- find_coordinated_programs(fit, n_perm = 200, seed = 1)
plot_cmp_nodes(cmp)
```
