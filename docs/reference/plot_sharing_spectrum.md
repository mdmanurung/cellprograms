# Plot the sharing spectrum

Per-pair summary of \`sharing_spectrum\`: the smallest principal angle
between each cell-type pair's program-activity subspaces, point size =
number of shared subspace dimensions (cosine \> 0.5), color =
permutation test outcome. Hollow points are underpowered pairs.

## Usage

``` r
plot_sharing_spectrum(spectrum)
```

## Arguments

- spectrum:

  Output of \`sharing_spectrum\`.

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
plot_sharing_spectrum(sharing_spectrum(fit, n_perm = 200))
```
