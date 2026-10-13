# Plot artifact-gene flags

Bar plot of the fraction of each program's top genes that are V(D)J or
donor-genotype-marker genes, with the flag threshold marked.

## Usage

``` r
plot_artifact_flags(flags, threshold = 0.2)
```

## Arguments

- flags:

  Output of \`flag_artifact_genes\`.

- threshold:

  Flag threshold (drawn as a dashed line; default 0.2).

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
plot_artifact_flags(flag_artifact_genes(fit))
```
