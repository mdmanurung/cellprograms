<div class="page-header">

# Plot program sharing results

<div class="d-none name">

`plot_sharing.Rd`

</div>

</div>

<div class="ref-description section level2">

Heatmap of pairwise cross-cell-type program correlations (the correlation pre-screen in \`program_sharing\`), with each program's final mashr-based classification appended to the axis labels.

</div>

<div class="section level2">

## Usage<a href="#ref-usage" class="anchor" aria-label="anchor"></a>

``` r
plot_sharing(share, level = c("program", "cell_type"), q_alpha = 0.05)
```

</div>

<div class="section level2">

## Arguments<a href="#arguments" class="anchor" aria-label="anchor"></a>

share<a href="#arg-share" class="anchor" aria-label="anchor"></a><br>
Output of \`program_sharing\`.

level<a href="#arg-level" class="anchor" aria-label="anchor"></a><br>
\`"program"\` (per-program matrix; readable up to a few dozen programs) or \`"cell_type"\` (mean pre-screen correlation per cell-type pair; the cohort-scale view).

q_alpha<a href="#arg-q-alpha" class="anchor" aria-label="anchor"></a><br>
Significance level for marking a match (default 0.05).

</div>

<div class="section level2">

## Value<a href="#value" class="anchor" aria-label="anchor"></a>

A ggplot object.

</div>

<div class="section level2">

## Examples<a href="#ref-examples" class="anchor" aria-label="anchor"></a>

``` r
sim <- simulate_pseudobulk(n_samples = 24, n_genes = 300, program_size = 30, seed = 1)
x <- as_cell_program_data(sim$pseudobulk)
fit <- canonicalize_programs(fit_celltype_programs(x, max_factors = 3, seed = 1))
if (requireNamespace("mashr", quietly = TRUE)) { # program_sharing needs mashr (Suggests)
  share <- program_sharing(fit, n_perm = 200, seed = 1)
  plot_sharing(share)
}
```

</div>
