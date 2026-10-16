# Plot program sharing results

Heatmap of pairwise cross-cell-type program correlations (the
correlation pre-screen in \`program_sharing\`), with each program's
final mashr-based classification appended to the axis labels.

## Usage

``` r
plot_sharing(share, level = c("program", "cell_type"), q_alpha = 0.05)
```

## Arguments

- share:

  Output of \`program_sharing\`.

- level:

  \`"program"\` (per-program matrix; readable up to a few dozen
  programs) or \`"cell_type"\` (mean pre-screen correlation per
  cell-type pair; the cohort-scale view).

- q_alpha:

  Significance level for marking a match (default 0.05).

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
if (requireNamespace("mashr", quietly = TRUE)) { # program_sharing needs mashr (Suggests)
  share <- program_sharing(fit, n_perm = 200, seed = 1)
  plot_sharing(share)
}
```
